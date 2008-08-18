#include <iostream>
#include <cmath>

#include "RKIntegrators.h"
#include <Numeric/arrayobject.h>
#define PyArray_DATA(obj) ((void *)(((PyArrayObject *)(obj))->data))
#include "triple.h"
#include <boost/python/extract.hpp>
#include <boost/python/to_python_converter.hpp>

const double RKIntegrator::speedOfLight = 2.99792458e8; // m/sec
const double RKIntegrator::speedOfLightSq = RKIntegrator::speedOfLight * 
                                           RKIntegrator::speedOfLight; // (m/sec)^2
const double RKIntegrator::electronCharge = -1.602176e-19; // Coulomb
const double RKIntegrator::electronMass = 9.1093826e-31;
const double RKIntegrator::electronMassEV = 510998.92;
const double RKIntegrator::ChargeOverMass = RKIntegrator::electronCharge/
                                            RKIntegrator::electronMass;
// To be studied a bit: which RK to use... 
const gsl_odeiv_step_type * RKIntegrator::myOdeiv_step_type = gsl_odeiv_step_rk8pd;

static const double RKIntegratorECloudCCChargeOverMass = -1.75882e11;

int RKIntegrator::signChange=1;

RKIntegrator::RKIntegrator(bool isRelativistic):
isRel(isRelativistic),
trajToFile(false),
trajFName("RKIntergatorTrajectoryDump.txt"),
precisionStep(1.0e-6),
stepRatio(5.),
inBeamPipe(true),
theDebugFlag(true),
errorInStep(false),
gamProtonBunch(8.938/0.938),
gamProtonBunchSq(gamProtonBunch*gamProtonBunch),
staticFieldModel(UNIFORM),
BFieldStrength(5.e-4),
BFieldMaxNonUnif(5.e-4),
BFieldMinNonUnif(5.e-4),
theTAbsolute(0.),
theTStart(0.),
theDeltaTGoal(0.),
theCurrentPotential(0.),
relDynamic(true),
maxXBeamPipe(5.e10),
maxYBeamPipe(5.e10),
minMaxRPipe(std::min(maxXBeamPipe, maxYBeamPipe)),
nRecursive(0),
relEFieldChange(0.01),
potentialUnits(1.0),
eFieldUnits(1.0),
thefOutPtr(0),
MIMagnetData(0)
{
  
  sOdeivStepAlloc = gsl_odeiv_step_alloc (myOdeiv_step_type, 6);
  cOdeivControl = gsl_odeiv_control_y_new (precisionStep, 0.0); 
  eOdeivEvolve = gsl_odeiv_evolve_alloc (6); 
  for (int k=0; k != 3; k++) BFieldStatic[k]=0.;
  for (int k=0; k != 6; k++) theVect6D[k]=0.;
//  std::cerr << " RKIntegrator instantiated!... " << std::endl;
}

RKIntegrator::~RKIntegrator() {
  gsl_odeiv_evolve_free (eOdeivEvolve);
  gsl_odeiv_control_free (cOdeivControl);
  gsl_odeiv_step_free (sOdeivStepAlloc);
  if (trajToFile && (thefOutPtr != 0)) this->closeTrajectoryFile();
  delete MIMagnetData;
}

void RKIntegrator::reOpenTrajectoryFile(const char *fName) 
{ 
    if (thefOutPtr != 0) this->closeTrajectoryFile();
    if (fName == NULL) trajToFile = false; 
    else {
       trajFName=std::string(fName); trajToFile=true; 
       thefOutPtr = new std::ofstream(trajFName.c_str());
       writeFOutTrajPreamble();
   } 
}

void RKIntegrator::reOpenTrajectoryFilePy(PyObject *fName) {
  char const* ff = boost::python::extract<char const*>(fName);
  this->reOpenTrajectoryFile(ff);
}
void RKIntegrator::setTrajectoryFileNamePy(PyObject *fName) {
  char const* ff = boost::python::extract<char const*>(fName);
  if (ff == 0) {
    std::cerr << " Null file name pointer, quit !! " << std::endl; exit(2);
  }
  this->setTrajectoryFileName(ff);
}
// Boost Magic.. 
PyObject *RKIntegrator::getTrajectoryFileNamePy() {
  return boost::python::incref(boost::python::object(trajFName).ptr());
}

void RKIntegrator::writeFOutTrajPreamble() {
  if (thefOutPtr == 0) return;
  // time, 6D phase space, with bx = beta_x.., Kinetic energy, Maxwell Field E, B, and potential 
  (*thefOutPtr) << " t x bx y by z bz eK Ex Ey Ez Bx By Bz VP" << std::endl; 

}	 

int RKIntegrator::propagateF(double *vect6D, const double *fields, double deltaT) {

 // Find out if we are in the Non-relativistc Regime (for electrons) 
 // Not yet implemented.. Done through the logical flag isRel 
  thefOutPtr=0;
  checkBeamPipeBoundary(vect6D);
  errorInStep=false;
  if (!inBeamPipe) return 0;
  nRecursive = 0;
  theTAbsolute = 0.; theDeltaTGoal = deltaT; theTStart = 0.;
  theCurrentPotential = 0.; // fixing the gauge, O.K., as the time is fixed... 
  if (trajToFile) { 
    thefOutPtr = new std::ofstream(trajFName.c_str());
    writeFOutTrajPreamble();
  }
 
 for (int k=0; k !=6; ++k) theVect6D[k] = vect6D[k];
 int nStep = this->propagateStepField(fields, deltaT);
 for (int k=0; k !=6; ++k) vect6D[k] = theVect6D[k];
 if (trajToFile) { 
   thefOutPtr->close();
   delete thefOutPtr;
   thefOutPtr=0;
 }
 return nStep;
}

int RKIntegrator::propagateFPy(PyObject *vect6DPy, PyObject *fieldsPy, double deltaT) {
    double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
    const double *fields = reinterpret_cast<double*>(PyArray_DATA(fields));
    this->propagateF(loc, fields, deltaT);
}
int RKIntegrator::propagateV(double *vect6D, const Real_scalar_field &potential, 
                              double tStart, double deltaT, double *tFinal) {

  bool debugFlag=theDebugFlag;
  errorInStep=false;
  if (debugFlag) 
     std::cerr << " Entering RKIntegrator::propagateV ... " << std::endl; 
  checkBeamPipeBoundary(vect6D);
  if (!inBeamPipe) return 0;
  // Already done! 
//  if ((thefOutPtr == 0) && trajToFile) { 
//    thefOutPtr = new std::ofstream(trajFName.c_str());
//    writeFOutTrajPreamble(); 
//  }
  // Crude adaptive stepping through a time dependend field 
  // Compute the field where we are...  (i) 
  // Transfor in lab frame, assuming electron velocity is small.. 
  // linear step, 1% of deltaT
  // Compute the field where we'll be, at t + 0.01*delta t 
  // if (relative value of field greater than 1%, ), divide the step... 
  // set the step, move. 
  // Go back to (i) if..
  // !! We currently set the Bfield to 0. 
  double locationS[3];
  double locationN[3]; 
  double fieldHere[6], fieldThere[6];
  double ETmp1[6], ETmp2[6];
  double tEnd = tStart + 0.999*deltaT;
  double deltaTCurrentExpected = deltaT;
//  std::cerr << " deltaTCurrentExpected = " << deltaTCurrentExpected << std::endl;
  double tNow = tStart;
  double epsilField = 1.0e10;
  const double epsilFieldTarget = relEFieldChange;
  std::vector<double> gSize=potential.get_physical_size();
  // Set the minimum time step in the loop below. 
  // If the grid is define too coarsely, not satisfactory convergence 
  // on local homogeneity can be set... Set to 1% of the bunch length.. 
  double minTimeStep=0.001*gSize[2]/speedOfLight;
  int nStep=0;
  theTAbsolute = tStart; theDeltaTGoal = deltaT;
  theTStart = tStart;
  for (int k=0; k!=6; k++) theVect6D[k] = vect6D[k]; 
  while (tNow < tEnd) {
    double deltaTCurrent = deltaTCurrentExpected;
    epsilField = 1.0e10;
    while (epsilField > epsilFieldTarget) {
      if (debugFlag) { 
        std::cerr << " Checking gradient at ";
        for (int k=0; k !=6; k++) std::cerr << ", " << theVect6D[k]; std::cerr << std::endl; 
        std::cerr << " .... deltaTCurrent " << deltaTCurrent << std::endl;
      }  
      for (int k=0; k !=3; k++) locationS[k] = theVect6D[2*k];
      // Beam moves along Z, neglect transverse motion over the length of one bunch..
      locationS[2] += speedOfLight*tNow;   
      // Compute the difference of potential.  If too small, field assumed to vanish 
      double potHere = potentialUnits*potential.get_val(locationS);
      for (int k=0; k !=3; k++) {
        ETmp1[k] = -1.*eFieldUnits*potential.get_deriv(locationS, k);
        ETmp1[k+3] = 0.; // No magnetic field from the proton bunch..., 
	                 // in the reference of the bunch
      }
      updateBField();  
      this->fieldBunchTransBeamToLab(ETmp1, fieldHere);
      for (int k=0; k !=3; k++) 
          locationN[k] = theVect6D[2*k] + theVect6D[2*k+1]*speedOfLight*deltaTCurrent;  
      locationN[2] += speedOfLight*(tNow+deltaTCurrent);   
      for (int k=0; k !=3; k++) {
        ETmp2[k] = -1.0*eFieldUnits*potential.get_deriv(locationN, k);
        ETmp2[k+3] = 0.;
      }
      double potThere = potentialUnits*potential.get_val(locationN);
      // Protect against indefinite fields (zero potentials) 
      if ((std::abs(potHere) < 0.01) && (std::abs(potThere) < 0.01)) {
         for (int k=0; k !=6; k++) {
	   fieldHere[k]=0.;
//	   if (k > 2) fieldHere[k] += BFieldStatic[k-3];
	   if (k > 2) fieldHere[k] += BFieldFromModel[k-3];
	 }
         break;
      }
      double deltaPot = std::abs(potHere - potThere);
      this->fieldBunchTransBeamToLab(ETmp2, fieldThere);
      double rDiffMax = 0.; double relDiff[3];
      for (int k=0; k !=3; k++) {
         if (std::abs(ETmp1[k]) > 1.0) { // below 1 Volt, we don't care... 
           relDiff[k] = std::abs((fieldThere[k] - fieldHere[k])/fieldHere[k]);
	   if (relDiff[k] > rDiffMax) rDiffMax = relDiff[k]; 
         }
      }
      if (debugFlag) {
        std::cerr << " relDiff ";  
        for (int k=0; k !=3; k++) std::cerr << ", " << relDiff[k]; std::cerr << std::endl;
      }  
      epsilField = rDiffMax;
      // Check that we won't go relativistic in the integration... The energy can go above 
      // 0.5 the electron mass (To be tuned!).  This, for all three direction
       
      for (int k=0; k!=3; k++) {
        double deltaE = std::abs(locationS[k] - locationN[k])*deltaPot;
        if (deltaE > 0.5*electronMassEV) {
	   epsilField=1.0e10;
	   std::cerr << " Extraordinarly strong field " << std::endl;
	   std::cerr << " Delta Potential " << deltaPot << " Along axis " 
	          << k << " delta l " <<  std::abs(locationS[k] - locationN[k])  
		  << " deltaE " << deltaE << std::endl;
	 }
      } 
      
      if (epsilField < epsilFieldTarget)  break;
//      std::cerr << " ..Relative difference " 
//                << epsilField << " deltaTCurrent " << deltaTCurrent << std::endl;
      deltaTCurrent /= 2; 
      if (deltaTCurrent < minTimeStep) break; // Potential just not smooth enough, gave up 
    } // on adjusting the step.. 
    // Now intergate over this step..
    theDeltaTGoal = deltaT;
    nRecursive = 0;
    double potS = potentialUnits*potential.get_val(locationS);
    theCurrentPotential = potS; 
    if (debugFlag) { // debugging... 
      double potN = potentialUnits*potential.get_val(locationN); 
      std::cerr << "---------------" << std::endl << " Decided on a step, start, vect6D ";
      for (int k=0; k !=6; ++k) std::cerr << ",  " <<  theVect6D[k];
      std::cerr << std::endl << " LocationS ";
      for (int k=0; k !=6; ++k) std::cerr << ",  " <<  locationS[k];
       std::cerr << std::endl << " .... Potential Here, Proton Bunch ref = "
	 << potS << " there " 
	 << potN << " diff " <<  (potS - potN) << std::endl;
       std::cerr << " .... fieldHere ";
      for (int k=0; k !=6; ++k) std::cerr << ", " <<  fieldHere[k];
      std::cerr << std::endl << " .... fieldThere  ";
      for (int k=0; k !=6; ++k) std::cerr << ", " <<  fieldThere[k];
      std::cerr << std::endl << " deltaTCurrent " << deltaTCurrent 
                << " tNow " << tNow << std::endl;
    }
    for (int k=0; k !=6; ++k) theVect6D[k] = vect6D[k];
    nStep += this->propagateStepField(fieldHere, deltaTCurrent);
    tNow = theTAbsolute;
    for (int k=0; k !=6; ++k) vect6D[k] = theVect6D[k];
//    std::cerr << " After step num " << nStep << ", end, vect6D ";
//    for (int k=0; k !=6; ++k) std::cerr << ",  " <<  vect6D[k];
//    std::cerr << std::endl;
    if (!inBeamPipe) break;
    // Try to upgrade the step size .. 
    deltaTCurrentExpected = 2*deltaTCurrent;
//    if (nStep > 20) {
//      std::cerr << " And quit after few steps in RKIntegrator::propagateV " << std::endl;
//      std::cerr << " tNow " << tNow << " tEnd " << tEnd << std::endl;
//      exit(2);
//    }
//    
 } // on step driven by field uniformity. 
 *tFinal = tNow; 
 return nStep;
 
}

int RKIntegrator::propagateVPy(PyObject *vect6DPy, Real_scalar_field &phi, 
                 double tStart, double deltaT, PyObject *timePy) {
   double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
   double *tFinal = reinterpret_cast<double*>(PyArray_DATA(timePy));		 
   this->propagateV(loc, phi, tStart, deltaT, tFinal);		 
}

int RKIntegrator::propThroughBunch(double *vect6D,  const Real_scalar_field &potential, 
                 double tOffset, double *tF){

  errorInStep=false;
  double tFTmp = tOffset;
  double tStart = tOffset;
  double rIn = std::sqrt(vect6D[0]*vect6D[0] + vect6D[2]*vect6D[2]); 
//  std::cerr << " RKIntegrator::propThroughBunch, rIn, C++ " << rIn << std::endl;
  int nStep = 0;
  checkBeamPipeBoundary(vect6D);
  if (!inBeamPipe) {
//    std::cerr <<  " PropThroughBunch: electron already out of Beam Pipe " << std::endl;
    return 0;
  }
  std::vector<double> gSize=potential.get_physical_size();
  // Assume the largest size is in Z... 
//  double mxSize=-999999.;
//  for (int k=0; k != gSize.size(); ++k)
//    if (std::abs(gSize[k]) > mxSize) mxSize =  gSize[k];
//  Bunch alway longer than wide.. 
//  mxSize =  gSize[2];
//  std::cerr << " RKIntegrator::propThroughProtonBunch,sizes, x = " << gSize[0] << 
//              ",y =  " << gSize[1] << ", z = " << gSize[2] << std::endl;
//  std::cerr << " And quit " << std::endl; exit(2);
  double betaBunch = std::sqrt(1. - 1./gamProtonBunchSq);
  // Increased a bit, by 0.5%, divide by two..
//  mxSize *= 0.5025;
  // Subtract from where we currently are.. No! done in the upstream class...  
//  double dz = mxSize - vect6D[4];
  // Assume the largest size is in Z... 
  double tLast = gSize[2]/(betaBunch*speedOfLight);
//  std::cerr << " RKIntegrator::propThroughProtonBunch, mxSize = " << mxSize
//            << " zinit " << vect6D[4] << " tLast " << tLast << std::endl;
  int kStep = 0;	 
  while (tFTmp < tLast) {
   if (theDebugFlag)    std::cerr << " propThroughBunch: Macro-step number kStep " 
              << kStep << std::endl;
    nStep += this->propagateV(vect6D,  potential, tStart, 5.0e-10, &tFTmp);
    if (!inBeamPipe) break;
    tStart = this->getTime();
    kStep++;
  }
  (*tF) = tFTmp;
  return nStep;		 
}

int RKIntegrator::propThroughBunchPy(PyObject *vect6DPy, Real_scalar_field &phi, 
                 double tOffset, PyObject *tFPy) {
  		 
  double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
//  double rIn = std::sqrt(loc[0]*loc[0] + loc[2]*loc[2]); 
//  std::cerr << " RKIntegrator::propThroughBunch, rIn, C++ " << rIn << std::endl;   
   double *tFinal = reinterpret_cast<double*>(PyArray_DATA(tFPy));
   int nn = this->propThroughBunch(loc, phi, tOffset, tFinal);
   return nn;		 
}

int RKIntegrator::propBetweenBunches(double *vect6D, 
                         double tOffset, double maxTime, double *tFinal) {

  errorInStep=false;
  checkBeamPipeBoundary(vect6D);
  if (!inBeamPipe) return 0;
  // Norm of the field 
  double bNormSq = 0.;
  double betaSq = 0.; 
  for (int k=0; k!=3; k++) {
    bNormSq += BFieldMaxNonUnif*BFieldMaxNonUnif;
    betaSq += vect6D[2*k+1]*vect6D[2*k+1];
  }
  double bNorm = std::sqrt(bNormSq);
  double beta = std::sqrt(betaSq);
  if (theDebugFlag)
    std::cerr << " propBetweenBunches: bNorm " << bNorm << " beta " << beta << std::endl;
  double gamma = std::sqrt(1.0/(1.0-betaSq));
  double e=electronMassEV*gamma*1.0e-9; // need energy in GeV to compute rho.. 
  double p=beta*e; // in GeV/c			 
  // Compute the Larmor radius			 
  double cost=0.;
  updateBField();
  for (int k=0; k !=3; k++) cost += vect6D[2*k+1]*BFieldFromModel[k];
  cost /= (bNorm*beta); 
  double pEff = p*std::sqrt(1.0 - cost*cost);
  if (bNorm < 0.1) bNorm = 0.5e-4; // minimum field, to avoid unphysical Nan... 
  double rho = pEff/(0.3*bNorm);
  double dt = rho/(beta*speedOfLight); // one steradian turn..
  double tau = 2.0*M_PI*rho/(beta*speedOfLight);
  double vParallel = speedOfLight*beta*cost;
  double timeToWall = minMaxRPipe/vParallel; // maximum time to wall.. 
  if (theDebugFlag) std::cerr << " .... Larmor radius " 
            << rho << " dt " << dt << " Freq. " 
	    << (1.0/tau) << " App. time to Wall " << timeToWall <<  std::endl;
	    
  theTAbsolute=tOffset;
  int nStep = 0;
  int nCallProp=0;
  double fieldHere[6];
  for (int k=0; k !=3; ++k) fieldHere[k] = 0.; 
//  for (int k=0; k !=3; ++k) fieldHere[k+3] = BFieldStatic[k];
  for (int k=0; k !=3; ++k) fieldHere[k+3] = BFieldFromModel[k];
  for (int k=0; k != 6; ++k) theVect6D[k] = vect6D[k]; 
  double vect6DPrev[6]; 
  double dz = 1000.0*precisionStep;
  // Set the precision to boundary to 10 time the precision of 
  while ((theTAbsolute < (tOffset+maxTime)) && (dz > 10.0*precisionStep)) {
    for (int k=0; k != 6; ++k) vect6DPrev[k] = theVect6D[k];
    if (theDebugFlag) {
      std::cerr << " dt is now " << dt << " tAbs " 
             << theTAbsolute << std::endl << " vect6D " ;
       for (int k=0; k !=6; k++) std::cerr << ", " << theVect6D[k]; std::cerr << std::endl;
    }
    theTStart = theTAbsolute; // by pass a previously installed chonology check. 
    nRecursive=0;
    updateBField();
    for (int k=0; k !=3; ++k) fieldHere[k+3] = BFieldFromModel[k];
    nStep += propagateStepField(fieldHere, dt);
    nCallProp++;
    // Geometry
    if (!inBeamPipe) { // Go back a few steps... Brute force bisecting.. 
      theTAbsolute -= dt;
      for (int k=0; k != 6; ++k) theVect6D[k] = vect6DPrev[k];
      dt /=2.;
      dz = dt*beta*speedOfLight;    
    }
    if (theTAbsolute < tOffset) break; // something wrong !
    if (nCallProp > 10000) {
      std::cerr << " prop between bunches Greater than 10000 integration steps.. " << std::endl;
      theTAbsolute = tOffset+maxTime+0.7e-12; // jump in time coordinate, we are going nowhere..
      break; 
    }  
  }
 *tFinal = theTAbsolute;
 for (int k=0; k != 6; ++k) vect6D[k] = theVect6D[k];
// std::cerr << " Quit from prop between bunches, after first track... " << std::endl;
// exit(2);
 return nStep; 
}

int RKIntegrator::propBetweenBunchesPy(PyObject *vect6DPy, 
                  double tOffset, double maxTime, PyObject *tFPy) {
		  
   double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
   double *tFinal = reinterpret_cast<double*>(PyArray_DATA(tFPy));		 
   return this->propBetweenBunches(loc, tOffset, maxTime, tFinal);		  
}

void RKIntegrator::closeTrajectoryFile() {
 if (trajToFile && (thefOutPtr != 0)) { 
   thefOutPtr->close();
   delete thefOutPtr;
   thefOutPtr = 0;
   trajFName=std::string("");
   trajToFile=false;
 }
}

int RKIntegrator::propagateStepField(const double *fields, double deltaT) {

  if (theDebugFlag) {
    std::cerr << " Entering propagateStepField, request time step " 
            << deltaT << " Current time " << theTAbsolute 
	    << " beta Z " << theVect6D[5] << std::endl;
    std::cerr << " Ey " << fields[1] << " nRecursive " << nRecursive  <<  std::endl; 
  }
  if (nRecursive > 10) {
    std::cerr << " Excessive recursion in RKIntegrator::propagateStepField, stop here !" << std::endl;
    errorInStep=true;
    return 0;
  }
  if (!checkBeamPipeBoundary(theVect6D)) return 0; 
  // Check the end game for recursion... 
  if ((theTAbsolute-theTStart) > theDeltaTGoal) return 0; // And no action!. 
  double t=theTAbsolute;
  // Extra copy because of recusrsion... 
  // (and automatic switching Relativistic not non-rel. )
  double vect6D[6]; for (int k=0; k !=6; k++) vect6D[k] = theVect6D[k];
//  std::cerr << " .... start, vect6D ";
//    for (int k=0; k !=6; ++k) std::cerr << ",  " <<  vect6D[k];
//  std::cerr << std::endl;
  gsl_odeiv_system sys;
  if (!isRel) {
    sys.function = RKIntegrator::funcNR; 
    sys.jacobian=  RKIntegrator::jacNR;
  } else { 
    sys.function =  RKIntegrator::funcR;
    sys.jacobian =  RKIntegrator::jacR;
//    std::cerr << " Going relativistic.. " << std::endl;
  }
  sys.dimension = 6;
  sys.params = (double *) fields;
  
  double t1 = theTAbsolute+deltaT; double h = precisionStep;
  double dt1 = deltaT/stepRatio; // arbtrary for now.... 
  // Go to velocities in m/Sec (SI units in this integrator) 
  for (int k=1; k!= 7; k+=2) vect6D[k] *= speedOfLight;
  // Iterate. 
  int n = 0;
//  std::cerr << " t = " << t << " t1 " << t1 << std::endl;
  
  double tBefore = theTAbsolute;
  bool goingForward = true; // Sometimes, the integrator make a step back....
  while (t < t1)
    {
      if (trajToFile) {
        double betax = vect6D[1]/speedOfLight; double betay = vect6D[3]/speedOfLight; 
        double betaz = vect6D[5]/speedOfLight; 
        (*thefOutPtr) << " " << (theTAbsolute)*1.0e9 << " " << vect6D[0]*1.0e3 << " " << betax;
        (*thefOutPtr) << " " << vect6D[2]*1.0e3 << " " << betay;
        (*thefOutPtr) << " " << vect6D[4]*1.0e3 <<" " << betaz;
        double betaSq = betax*betax + betay*betay + betaz*betaz;
        double eK = electronMassEV*betaSq/2.;
        if (betaSq > 0.0001) {
          double gam=std::sqrt(1.0/(1.0 - betaSq));
	  eK = electronMassEV*(gam - 1.0);
        }
        (*thefOutPtr) << " " << eK;  
	for (int k=0; k!=6; ++k) (*thefOutPtr) << " " << fields[k];
	(*thefOutPtr) << " " << theCurrentPotential;
        (*thefOutPtr) << std::endl;  
      }
      // Set the required precision... 
     {
         double betax = vect6D[1]/speedOfLight; double betay = vect6D[3]/speedOfLight; 
         double betaz = vect6D[5]/speedOfLight; 
         double betaSq = betax*betax + betay*betay + betaz*betaz;
         double gammaSq = 1.0/(1.-betaSq);
	 double fact = betaSq*gammaSq;
	 if (fact > 1.0e10) {
	   std::cerr << " Setting new precision.., betagam sq " << fact << std::endl;
           double aPrec = precisionStep/(fact*fact);
	   this->setPrecisionStep(aPrec);
	 }
// Check if the equation of motion are appropriate... Relativistic vs non-relativistic... 
//

         if (isRel && relDynamic) {
           if (betaSq < 1.0e-4) { // we can go non-relativistic... 
             // Revert to beta units (SI units in this integrator) 
             for (int k=1; k!= 7; k+=2) vect6D[k] /= speedOfLight;  
             gsl_odeiv_evolve_reset(eOdeivEvolve);  
	     double deltaTNext = deltaT-(theTAbsolute-theTStart);
	     isRel = false;
//	     std::cerr << " Going non-relativistic...at time  " << theTAbsolute << std::endl;
//	     std::cerr << " x =  " << vect6D[0] << " y " 
//	               << vect6D[2] << " z " << vect6D[4] << " betaSq " << betaSq << std::endl;
             for (int k=0; k !=6; k++) theVect6D[k] = vect6D[k];
	     nRecursive++;		       
	     this->propagateStepField(fields, deltaTNext);
	     return n;
//	      
//             if ((theTAbsolute-theTStart) > theDeltaTGoal) return 0;
//	     if (!inBeamPipe) return 0;
	        // No action, we done... Recursion end game. 
	   }
         } else if (!isRel && relDynamic) { // we are non-relativistic.. Place a tight limit, accuracy 
	          // matters more than speed.  
           if (betaSq > 1.0e-4) { // we must go relativistic... 
             // Revert to beta units (SI units in this integrator) 
             for (int k=1; k!= 7; k+=2) vect6D[k] /= speedOfLight;  
             gsl_odeiv_evolve_reset(eOdeivEvolve);  
	     double deltaTNext = deltaT-(theTAbsolute-theTStart);
	     isRel = true;
//	     std::cerr << " Going Relativistic at time " << theTAbsolute << std::endl;
//	     std::cerr << " x =  " << vect6D[0] << " y " 
//	               << vect6D[2] << " z " << vect6D[4] << " betaSq " << betaSq << std::endl;
             for (int k=0; k !=6; k++) theVect6D[k] = vect6D[k];		       
	     nRecursive++;		       
	     this->propagateStepField(fields, deltaTNext); 
	     return n;
//             if ((theTAbsolute-theTStart) > theDeltaTGoal) return 0;
//	     if (!inBeamPipe) return 0;
	   }
	 }
//         if (relDynamic & (n < 20)) 
//	       std::cerr << " at t Abs " << theTAbsolute << " betaSq " << betaSq << std::endl;
      }
      double tGoal = t + dt1; // for this sub-step... 
      tBefore = t;
//      if (isRel) std::cerr << " Rel.."; else std::cerr << " NRel..";
//      std::cerr << " Bf Step ... ";
//      for (int k=0; k !=6; ++k) std::cerr << ",  " <<  vect6D[k];
//      std::cerr << std::endl;
      int status = gsl_odeiv_evolve_apply (eOdeivEvolve, cOdeivControl, sOdeivStepAlloc,
                                           &sys, 
                                           &t, tGoal,
                                           &h, vect6D);
//      std::cerr << " t = " << t << std::endl;
      theTAbsolute += t - tBefore;
//      if (isRel) std::cerr << " Rel.."; else std::cerr << " NRel..";
//      std::cerr << " dt1 ... " << dt1 << " delta Reached " 
//                << t - tBefore << " tAbsolute= " << theTAbsolute << " vY " << vect6D[3] << std::endl;
//      std::cerr << " Af Step ... ";
//      for (int k=0; k !=6; ++k) std::cerr << ",  " <<  vect6D[k];
//      std::cerr << std::endl;
      if (tBefore > t) {
        std::cerr << " Backward step!.... at t = " << tBefore 
	          << " fatal confusion for now... "  << std::endl;
	goingForward = false;
	exit(2); // nned fixing, if need be....
      }

      if (status != GSL_SUCCESS) {
        std::cerr << " RKIntegrator failed in propagateNR, at time " << t << std::endl;
	std::cerr << " state vector ";
	for (int k=0; k !=6; k++) std::cerr << " " << vect6D[k]; std::cerr << std::endl;
        break;
      }
      // Geometry 
      if (!checkBeamPipeBoundary(vect6D)) {
//         std::cerr << " Stopping propagation, on beamPipe..." << std::endl; 
         break;
      }
      n++;
      if (theDebugFlag) { 
         std::cerr << " PropagateStepField, at step " << n << std::endl;
	 std::cerr << "   vect6D "; for (int k=0; k!=6; k++) std::cerr << ", " << vect6D[k];
	 std::cerr << std::endl;
      }
      // Against infinite loops! Needs to report error. 
      if (n > 1000) {
//        std::cerr << " PropagateStepField: Too many interation, and far from the beam bipe.." <<
//	       std::endl;
	double dBP = distToBeamPipe(vect6D);
//        std::cerr << " Distance to beam bipe.." <<  dBP << std::endl;
	double dBPInit = distToBeamPipe(theVect6D);
//        std::cerr << " Distance to beam bipe, start " <<  dBPInit << std::endl;
        if (dBP < .05) break; // 1% spatial error possible reaching the beam pipe. 
//	 std::cerr << "   vect6D "; for (int k=0; k!=6; k++) std::cerr << ", " << vect6D[k];
//	 std::cerr << std::endl;
         errorInStep = true;
         break;
      }
//      if (n > 5) {
//        std::cerr << " Enough debugging, stop here " << std::endl; exit(2);
//      }
    }
  // Revert to beta units (SI units in this integrator) 
  for (int k=1; k!= 7; k+=2) vect6D[k] /= speedOfLight;  
  for (int k=0; k !=6; k++) theVect6D[k] = vect6D[k];		       
//  std::cerr << " ....nRecursive= " << nRecursive << "  end, vect6D ";
//    for (int k=0; k !=6; ++k) std::cerr << ",  " <<  vect6D[k];
//  std::cerr << std::endl;
  gsl_odeiv_evolve_reset(eOdeivEvolve); 
//  std::cerr << " Exiting PropagateFieldStep at t = " <<  theTAbsolute << std::endl;
  return n;

}

int RKIntegrator::funcNR (double t, const double y[], double f[], void *params) {

  double *Field = (double *)params;
  //0-2 is Ex, Ey, Ez, 3-5 is Bx, By, Bz
  // y[0-5] is x 'x y y' z z', units are SI, meter and meter/sec 
  f[0] = y[1];
  f[1] = RKIntegrator::ChargeOverMass*signChange*(Field[0] + y[3]*Field[5] - y[5]*Field[4]);
  f[2] = y[3];
  f[3] = RKIntegrator::ChargeOverMass*signChange*(Field[1] + y[5]*Field[3] - y[1]*Field[5]);
  f[4] = y[5];
  f[5] = RKIntegrator::ChargeOverMass*signChange*(Field[2] + y[1]*Field[4] - y[3]*Field[3]);
  return GSL_SUCCESS;


}

int RKIntegrator::funcR (double t, const double yIn[], double f[], void *params) {
    
  double *Field = (double *)params;
  //0-2 is Ex, Ey, Ez, 3-5 is Bx, By, Bz
  // coord[0-5] is x ux y uy z uz

  double Ex, Ey, Ez, Bx, By, Bz;
  double x=yIn[0]; double vx=yIn[1]; // in m/sec, SI units... 
  double y=yIn[2]; double vy=yIn[3];
  double z=yIn[4]; double vz=yIn[5];

  double c, velsq, betasq, gamma;


  Ex = Field[0]; Ey = Field[1]; Ez = Field[2];

  Bx = Field[3]; By = Field[4]; Bz = Field[5];

//  velsq = pow(vx,2.) +  pow(vy,2.) +  pow(vz,2.);
  velsq = vx*vx +  vy*vy +  vz*vz;
  betasq = velsq/speedOfLightSq;
  if (betasq > 1.) {
    std::cerr << " Going tachionic in funcR, at x " << x << " y " << y << " z " << z << std::endl;
    std::cerr << " ... vx " << vx << " vy " << vy << " vz " << vz;
    std::cerr << " ... Ex " << Ex << " Ey " << Ey << " Ez " << Ez << std::endl;
    std::cerr << " Fatal, stop simulation.. " << std::endl;
    exit(2);
  }
  gamma = 1.0/std::sqrt((1.-betasq));

  f[0] = vx;
  f[2] = vy;
  f[4] = vz;
  
  double gamma2 = gamma*gamma;
  double gam2 = gamma2;
  double m = electronMass;
  double vx2 = vx*vx;
  double vy2 = vy*vy;
  double vz2 = vz*vz;
  double qOverM = signChange*electronCharge/electronMass;
  double coeffx = qOverM/ ((1+gam2*vx2/speedOfLightSq)*gamma);
  f[1] = coeffx*(Ex+(vy*Bz-vz*By));
  double coeffy = qOverM/ ((1+gam2*vy2/speedOfLightSq)*gamma);
  f[3] = coeffy*(Ey+ (vz*Bx-vx*Bz));
  double coeffz = qOverM/ ((1+gam2*vz2/speedOfLightSq)*gamma);
  f[5] = coeffz*(Ez+ (vx*By-vy*Bx));
  
  return GSL_SUCCESS;
    
}   

int RKIntegrator::jacNR (double t, const double y[], double *dfdy, 
                double dfdt[], void *params) {
		
  double *Field = (double *)params;
  //0-2 is Ex, Ey, Ez, 3-5 is Bx, By, Bz
  // in this example q/m, c are= 1, and the equ of motion non-relativistic 
  // y[0-5] is x ux y uy z uz

  std::cerr << " In jacNR, t = " << t << std::endl;
  std::cerr << " Wrong stuff, forget it... " << std::endl;
  exit(2);

  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 6, 6);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 0, 4, 0.0);
  gsl_matrix_set (m, 0, 5, 0.0);

  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, 0.0);
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, RKIntegratorECloudCCChargeOverMass*signChange*Field[5]);
  gsl_matrix_set (m, 1, 4, 0.0);
  gsl_matrix_set (m, 1, 5, -RKIntegratorECloudCCChargeOverMass*signChange*Field[4]);

  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, 0.0);
  gsl_matrix_set (m, 2, 3, 1.0);
  gsl_matrix_set (m, 2, 4, 0.0);
  gsl_matrix_set (m, 2, 5, 0.0);

  gsl_matrix_set (m, 3, 0, 0.0);
  gsl_matrix_set (m, 3, 1, -RKIntegratorECloudCCChargeOverMass*signChange*Field[5]);
  gsl_matrix_set (m, 3, 2, 0.0);
  gsl_matrix_set (m, 3, 3, 0.0);
  gsl_matrix_set (m, 3, 4, 0.0);
  gsl_matrix_set (m, 3, 5, RKIntegratorECloudCCChargeOverMass*signChange*Field[3]);

  gsl_matrix_set (m, 4, 0, 0.0);
  gsl_matrix_set (m, 4, 1, 0.0);
  gsl_matrix_set (m, 4, 2, 0.0);
  gsl_matrix_set (m, 4, 3, 0.0);
  gsl_matrix_set (m, 4, 4, 0.0);
  gsl_matrix_set (m, 4, 5, 1.0);

  gsl_matrix_set (m, 5, 0, 0.0);
  gsl_matrix_set (m, 5, 1, RKIntegratorECloudCCChargeOverMass*signChange*Field[4]);
  gsl_matrix_set (m, 5, 2, 0.0);
  gsl_matrix_set (m, 5, 3, -RKIntegratorECloudCCChargeOverMass*signChange*Field[3]);
  gsl_matrix_set (m, 5, 4, 0.0);
  gsl_matrix_set (m, 5, 5, 0.0);

  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[3] = 0.0;
  dfdt[4] = 0.0;
  dfdt[5] = 0.0;
  return GSL_SUCCESS;
		
}		

int RKIntegrator::jacR (double t, const double coord[], double *dfdy, 
                double dfdt[], void *params){
		
  double *Field = (double *)params;
  //0-2 is Ex, Ey, Ez, 3-5 is Bx, By, Bz
  // in this example q/m, c are= 1, and the equ of motion non-relativistic 
  // y[0-5] is x ux y uy z uz

  double Ex, Ey, Ez, Bx, By, Bz;
  double x, vx, y, vy, z, vz;

  double c, m, q, velsq, betasq, gamma;
  
/*
  double dvxx, dvxy, dvxz;
  double dvyx, dvyy, dvyz;
  double dvzx, dvzy, dvzz;
*/
//  double Cf1, Cf2, Cf3, Cf4, Cf5, Cf6, Cf7, Cf8, Cf9, Cf10;
//  double Cf11, Cf12, Cf13, Cf14, Cf15, Cf16, Cf17, Cf18;


  Ex = Field[0]; Ey = Field[1]; Ez = Field[2];

  Bx = Field[3]; By = Field[4]; Bz = Field[5];

  x = coord[0]; vx = coord[1]; 
  y = coord[2]; vy = coord[3]; 
  z = coord[4]; vz = coord[5]; 

  c = speedOfLight;
  q = signChange*electronCharge;
  m = electronMass;

  velsq = vx*vx +  vy*vy +  vz*vz;
  betasq = velsq/speedOfLightSq;
  gamma = 1.0/std::sqrt((1.-betasq));
  std::cerr << "Relativistic in jacR .. betasq " << betasq <<" gamma "<< gamma << std::endl;

/* From Spentz prototype..   

  Cf1=pow(c,4); //To be improved, not need to invoke the power function, 
  Cf2=1/Cf1;
  Cf3=1/m;
  Cf4=pow(c,2);
  Cf5=pow(gamma,2);
  Cf6=pow(vx,3);
  Cf7=pow(gamma,4);
  Cf8=1/gamma;
  Cf9=pow(vy,2);
  Cf10=pow(vx,2);
  Cf11=pow(gamma,5);
  Cf12=Cf4*Ey*q-Bz*Cf4*q*vx;
  Cf13=pow(vy,3);
  Cf14=pow(gamma,3);
  Cf15=-Cf4*Ez*q;
  Cf16=-By*Cf4*q*vx;
  Cf17=-2*Cf4*m;
  Cf18=2*m*Cf10;

  dvxx = -Cf2*Cf3*(2*m*Cf6*Cf7+2*Cf4*m*vx*Cf5+(-By*Cf4*q*vx*vz+Bz*Cf4*q* 
					       vx*vy+Cf4*Ex*q*vx)*gamma);
  dvxy = -Cf2*Cf3*Cf8*(2*m*Cf10*vy* 
		       Cf11+(-By*Cf4*q*vy*vz+Bz*Cf4*q*Cf9+Cf4*Ex*q*vy)*Cf5-Bz*Cf1* 
		       q);
  dvxz = -Cf2*Cf3*Cf8*(2*m*Cf10*vz*Cf11+((Bz*Cf4*q*vy+Cf4*Ex*q)* 
					 vz+By*Cf4*q*Cf9+By*Cf4*q*Cf10-By*Cf1*q)*Cf5+2*By*Cf1* 
		       q);
  dvyx = -Cf2*Cf3*Cf8*(2*m*vx*Cf9*Cf11+(Bx*Cf4*q*vx*vz-Bz*Cf4*q* 
					Cf10+Cf4*Ey*q*vx)*Cf5+Bz*Cf1*q);
  dvyy = -Cf2*Cf3*(2*m*Cf13*Cf7+2*Cf4* 
		   m*vy*Cf5+(Bx*Cf4*q*vy*vz+Cf12*vy)*gamma);
  dvyz = -Cf2*Cf3*Cf8*(2*m* 
		       Cf9*vz*Cf11+(Cf12*vz-Bx*Cf4*q*Cf9-Bx*Cf4*q*Cf10+Bx*Cf1*q)*Cf5-2*Bx* 
		       Cf1*q);
  dvzx = Cf2*Cf3*Cf8*((2*m*vx*Cf9+2*m*Cf6-2*Cf4*m*vx)*Cf11+2* 
		      Cf4*m*vx*Cf14+(Bx*Cf4*q*vx*vy-By*Cf4*q*Cf10-Cf4*Ez*q*vx)*Cf5+By*Cf1* 
		      q);
  dvzy = Cf2*Cf3*Cf8*((2*m*Cf13+(Cf18+Cf17)*vy)*Cf11+2*Cf4*m*vy* 
		      Cf14+(Bx*Cf4*q*Cf9+(Cf16+Cf15)*vy)*Cf5-Bx*Cf1*q);
  dvzz = Cf2*Cf3*((2* 
		   m*Cf9+Cf18+Cf17)*vz*Cf7+(Bx*Cf4*q*vy+Cf16+Cf15)*vz*gamma);

*/

// From Maxima stringout command, repeat.... 
/*
[[dvxx = 0,
dvxy = -2*Bz*q/((2*bet2-1)*m),
dvxz = 2*By*q/((2*bet2-1)*m)],
[dvyx = 2*Bz*q/((2*bet2-1)*m),
dvyy = 0,
dvyz = -2*Bx*q/((2*bet2-1)*m)],
[dvzx = -2*By*q/((2*bet2-1)*m),
dvzy = 2*Bx*q/((2*bet2-1)*m),
dvzz = 0]];
// Clearly different!!!!! 
*/
// C++ implementation of above.. 
 double dvxx = 0.;
 double dvxy = -2*Bz*q/((2*betasq-1)*m);
 double dvxz = 2*By*q/((2*betasq-1)*m);
 double dvyx = 2*Bz*q/((2*betasq-1)*m);
 double dvyy = 0.;
 double dvyz = -2*Bx*q/((2*betasq-1)*m);
 double dvzx = -2*By*q/((2*betasq-1)*m);
 double dvzy = 2*Bx*q/((2*betasq-1)*m);
 double dvzz = 0.;

// std::cerr << " dvxy " << dvxy << std::endl;

  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 6, 6);
  gsl_matrix * M = &dfdy_mat.matrix; 
  gsl_matrix_set (M, 0, 0, 0.0);
  gsl_matrix_set (M, 0, 1, 1.0);
  gsl_matrix_set (M, 0, 2, 0.0);
  gsl_matrix_set (M, 0, 3, 0.0);
  gsl_matrix_set (M, 0, 4, 0.0);
  gsl_matrix_set (M, 0, 5, 0.0);

  gsl_matrix_set (M, 1, 0, 0.0);
  gsl_matrix_set (M, 1, 1, dvxx);
  gsl_matrix_set (M, 1, 2, 0.0);
  gsl_matrix_set (M, 1, 3, dvxy);
  gsl_matrix_set (M, 1, 4, 0.0);
  gsl_matrix_set (M, 1, 5, dvxz);

  gsl_matrix_set (M, 2, 0, 0.0);
  gsl_matrix_set (M, 2, 1, 0.0);
  gsl_matrix_set (M, 2, 2, 0.0);
  gsl_matrix_set (M, 2, 3, 1.0);
  gsl_matrix_set (M, 2, 4, 0.0);
  gsl_matrix_set (M, 2, 5, 0.0);

  gsl_matrix_set (M, 3, 0, 0.0);
  gsl_matrix_set (M, 3, 1, dvyx);
  gsl_matrix_set (M, 3, 2, 0.0);
  gsl_matrix_set (M, 3, 3, dvyy);
  gsl_matrix_set (M, 3, 4, 0.0);
  gsl_matrix_set (M, 3, 5, dvzz);

  gsl_matrix_set (M, 4, 0, 0.0);
  gsl_matrix_set (M, 4, 1, 0.0);
  gsl_matrix_set (M, 4, 2, 0.0);
  gsl_matrix_set (M, 4, 3, 0.0);
  gsl_matrix_set (M, 4, 4, 0.0);
  gsl_matrix_set (M, 4, 5, 1.0);

  gsl_matrix_set (M, 5, 0, 0.0);
  gsl_matrix_set (M, 5, 1, dvzx);
  gsl_matrix_set (M, 5, 2, 0.0);
  gsl_matrix_set (M, 5, 3, dvzy);
  gsl_matrix_set (M, 5, 4, 0.0);
  gsl_matrix_set (M, 5, 5, dvzz);

  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[3] = 0.0;
  dfdt[4] = 0.0;
  dfdt[5] = 0.0;
  return GSL_SUCCESS;
    
}   
//  assume pure electric field.. and vx=vy=0, vz only 
void RKIntegrator::fieldBunchTransBeamToLab(const double *fieldIn, 
                                            double *fieldOut) {
				    
 double betaZ = std::sqrt(1.0-1.0/gamProtonBunchSq);
 fieldOut[0] = gamProtonBunch*fieldIn[0]; fieldOut[1] = gamProtonBunch*fieldIn[1];
 fieldOut[2] = fieldIn[2];
 fieldOut[3] = BFieldFromModel[0] + 1.0*gamProtonBunch*(betaZ/speedOfLight)*fieldIn[1];
 fieldOut[4] = BFieldFromModel[1] + 1.0*gamProtonBunch*(betaZ/speedOfLight)*fieldIn[0];
 fieldOut[5] = BFieldFromModel[2];  

}

bool RKIntegrator::checkBeamPipeBoundary(const double *v) {

  inBeamPipe = true; // Benefit of the doubt!  
  if ((std::abs(v[0]) > maxXBeamPipe) || (std::abs(v[2]) > maxYBeamPipe)) {
      inBeamPipe = false;
      return false;
  }
  double rrWallSq = (v[0]*v[0]/(maxXBeamPipe*maxXBeamPipe)) + 
                   (v[2]*v[2]/(maxYBeamPipe*maxYBeamPipe));
  if (rrWallSq > 1.) {
      inBeamPipe = false;
      return false;
  }
  return true;

}

double RKIntegrator::distToBeamPipe(const double *v) {

  if ((std::abs(v[0]) > maxXBeamPipe) || (std::abs(v[2]) > maxYBeamPipe)) 
      return 0.;
  double rrWallSq = (v[0]*v[0]/(maxXBeamPipe*maxXBeamPipe)) + 
                   (v[2]*v[2]/(maxYBeamPipe*maxYBeamPipe));
  return std::sqrt(std::abs(1.0-rrWallSq)); //effective metric ! Relative..

}

void RKIntegrator::setUnits(double totalCharge, double units0) {

// Total charge x 1.0/(4.0*pi*eps0) (in SI units) 
// The factor .04954 is a mystery one. To be fixed later... 
//  potentialUnits = totalCharge*8.98755e9/0.04954;
  potentialUnits = totalCharge*8.98755e9*4.0*M_PI; // Check with ../StudyPotential 
  // Grid size must matter for the potential.. Guess work.. 
  eFieldUnits = potentialUnits/units0;
//  std::cerr << " Potential Units " << potentialUnits << " And quit " <<  std::endl;
//  exit(2); 
}

void RKIntegrator::setFieldModel(BFieldModel aModel, double strength) {
   
   staticFieldModel=aModel;
   // compute the minim and max value of the field over ~ 2.7 m ( 6 x bunch length) 
   switch (aModel) {
    case UNIFORM: {
        double bNorm = 0.;
        for (int k=0; k!=3; k++) bNorm += BFieldStatic[k]*BFieldStatic[k];
        BFieldMaxNonUnif=std::sqrt(bNorm); 
        BFieldMinNonUnif=BFieldMaxNonUnif; }
      break;
    
    case QUADRUPOLE: { // Uniform in z, not quad rotation..
      if (MIMagnetData == 0) MIMagnetData = new MIDipQuadEdge();
      MIMagnetData->setMIEnergy(strength);
      BFieldStrength= strength*MIMagnetData->BQuadStrOverMIEnergy;
      double bxMax=BFieldStrength*maxXBeamPipe;
      double byMax=BFieldStrength*maxYBeamPipe;
      BFieldMinNonUnif = 1.0e-3;  // at center, roughly. 
      BFieldMaxNonUnif = std::max(bxMax, byMax); 
      break; }
      
    case DIPOLEEDGE:
      if (MIMagnetData == 0) MIMagnetData = new MIDipQuadEdge();
      BFieldStrength= strength*MIMagnetData->BDipoleOverMIEnergy;
      // Strength is in fact the MI energy... in this context.. 
      MIMagnetData->setMIEnergy(strength);
      BFieldMaxNonUnif=strength; 
      BFieldMinNonUnif = 0.5e-4;  // Earth field, far away..  
      break;
      
    case QUADRUPOLEEDGE: {
        if (MIMagnetData == 0) MIMagnetData = new MIDipQuadEdge();
        MIMagnetData->setMIEnergy(strength);
        BFieldStrength= strength*MIMagnetData->BQuadStrOverMIEnergy;
        double bxMax=BFieldStrength*maxXBeamPipe;
        double byMax=BFieldStrength*maxYBeamPipe;
        BFieldMaxNonUnif=std::max(bxMax, byMax); 
        BFieldMinNonUnif = 0.5e-4;  
      }
      break;      
    case DIPQUADEDGES:
      {
        if (MIMagnetData == 0) MIMagnetData = new MIDipQuadEdge();
      // Strength is in fact the MI energy... in this context.. 
        MIMagnetData->setMIEnergy(strength);
      // For sake of consitency, 
        BFieldStrength= strength*0.011409; // MIDipQuadEdge class 
        BFieldMaxNonUnif=BFieldStrength; 
        BFieldMinNonUnif = 0.5e-4;  // Earth field, far away..
       }  
      break;
      
   }
} 

void RKIntegrator::setFieldModelPy(int aModel, double strength) {
  // No protection,, tired of computology.. 
  BFieldModel  am = (BFieldModel) aModel;
  this->setFieldModel(am, strength);
}

void RKIntegrator::updateBField() {
  
  switch (staticFieldModel) {
     case UNIFORM:
       for (int k=0; k!=3; k++) BFieldFromModel[k] = BFieldStatic[k]; // old model. 
       break;
     case QUADRUPOLE:
       BFieldFromModel[2] = 0.;
       BFieldFromModel[0] = BFieldStrength*theVect6D[0];
       BFieldFromModel[1] = BFieldStrength*theVect6D[2];
//       std::cerr << " At x= " << theVect6D[0] << " Bx " << BFieldFromModel[0] << 
//                " y = " << theVect6D[2] << " By " << BFieldFromModel[1] << std::endl;
       break;
    case DIPOLEEDGE:
      { 
        for (int k=0; k!=3; k++) BFieldFromModel[k] = 0.5e-4; // earth field. 
        double fact = MIMagnetData->getFieldYFact(theVect6D[4]);
        BFieldFromModel[1] = BFieldStrength*fact;
      }
      break;
      
    case QUADRUPOLEEDGE: 
      {
        double fact = MIMagnetData->getFieldYFact(theVect6D[4]);
        BFieldFromModel[1] = BFieldStrength*fact*theVect6D[2];
        BFieldFromModel[0] = BFieldStrength*fact*theVect6D[0];
        BFieldFromModel[2] = 0.5e-4;
      }
      break;
      
    case DIPQUADEDGES:
      { 
        double loc[3]; 
        for (int k=0; k!=3; k++) loc[k] = theVect6D[2*k];
        BFieldFromModel[2] = 0.5e-4;
        BFieldFromModel[0] = MIMagnetData->getFieldX(loc);
        BFieldFromModel[1] = MIMagnetData->getFieldY(loc);
      }
      break;
      
     default:
       for (int k=0; k!=3; k++) BFieldFromModel[k] = 1.0e-4;
       break;
    }

}
