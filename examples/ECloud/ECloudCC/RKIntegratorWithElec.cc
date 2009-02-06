#include <iostream>
#include <cmath>

#include "RKIntegrators.h"
#include <numpy/arrayobject.h>
#define PyArray_DATA(obj) ((void *)(((PyArrayObject *)(obj))->data))
#include "triple.h"
#include <boost/python/extract.hpp>
#include <boost/python/to_python_converter.hpp>

// Hit level routines to integrate the equation of motion of electrons 
// within the field of the proton bunch, and the self field of these electrons
// The calculation of the potential from these electron is done elsewhere, 
// i.e., we only deal here with one electron at a time. 
// 
// P.L. Aug 18 2008 

int RKIntegrator::propagateV(double *vect6D, const Real_scalar_field &potential,
                              const Real_scalar_field &electronPotential, 
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
  double locationSBP[3]; // not shifted, in BP.
  double locationNBP[3]; 
  double fieldHere[6], fieldThere[6];
  double ETmpHere[6], ETmpThere[6], ETmpHereElec[6], ETmpThereElec[6];
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
      for (int k=0; k !=3; k++) {
          locationS[k] = theVect6D[2*k];
	  locationSBP[k] = locationS[k]; 
      }	  
      // Beam moves along Z, neglect transverse motion over the length of one bunch..
      locationS[2] += speedOfLight*tNow;   
      // Compute the difference of potential.  If too small, field assumed to vanish 
      double potHere = potentialUnits*potential.get_val(locationS);
      double potHereElectron = potentialUnitsElectr*electronPotential.get_val(locationSBP);
      for (int k=0; k !=3; k++) {
        ETmpHere[k] = -1.*eFieldUnits*potential.get_deriv(locationS, k);
        ETmpHere[k+3] = 0.; // No magnetic field from the proton bunch..., 
	                 // in the reference of the bunch
      }
      updateBField();  
      this->fieldBunchTransBeamToLab(ETmpHere, fieldHere);
      // O.K., that was here.  Now check overt there, 
      // and see if it not too different.. 
      
      for (int k=0; k !=3; k++) {
          locationN[k] = theVect6D[2*k] + theVect6D[2*k+1]*speedOfLight*deltaTCurrent; 
	  locationNBP[k] =locationN[k];
      }
      locationN[2] += speedOfLight*(tNow+deltaTCurrent);   
      double potThere = potentialUnits*potential.get_val(locationN);
      double potThereElectron = potentialUnitsElectr*electronPotential.get_val(locationNBP);
      double deltaPot = std::abs(potHere - potThere);
      for (int k=0; k !=3; k++) {
        ETmpThere[k] = -1.0*eFieldUnits*potential.get_deriv(locationN, k);
        ETmpThere[k+3] = 0.;
      }
      this->fieldBunchTransBeamToLab(ETmpThere, fieldThere);
      
      // Before adding the field from electron, check if 
      // if the field from the proton bunch is O.K. 
      // Protect against indefinite fields (zero potentials) 
      if ((std::abs(potHere) < 0.01) && (std::abs(potThere) < 0.01)) {
         for (int k=0; k !=6; k++) {
	   fieldHere[k]=0.;
//	   if (k > 2) fieldHere[k] += BFieldStatic[k-3];
	   if (k > 2) fieldHere[k] += BFieldFromModel[k-3];
	 }
      }
      if ((std::abs(potHereElectron) > 0.01) || (std::abs(potThereElectron) > 0.01)) {
      // Sign to be checked!. 
        for (int k=0; k !=3; k++) {
          ETmpHereElec[k] = -1.0*eFieldUnitsElectr*electronPotential.get_deriv(locationSBP, k);
          ETmpHereElec[k+3] = 0.; // No magnetic field from these electrons  
          ETmpThereElec[k] = -1.0*eFieldUnitsElectr*electronPotential.get_deriv(locationNBP, k);
          ETmpThereElec[k+3] = 0.; // No magnetic field from these electrons  
        }      
      
        for (int k=0; k !=3; k++) {  // No magentic field from the electron, too low 
	                             // energy. 
	   fieldHere[k] += ETmpHereElec[k];
           fieldThere[k] += ETmpThereElec[k];
	}
      }
      
      double rDiffMax = 0.; double relDiff[3];
      for (int k=0; k !=3; k++) {
         if (std::abs(ETmpHere[k]) > 1.0) { // below 1 Volt, we don't care... 
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
      // Academic, for the most part.. And not even correct, since 
      // it is not estimated in the correct reference frame.. Leave as it...  
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
    double potSEl = potentialUnitsElectr*electronPotential.get_val(locationS); 
    theCurrentPotential = potS; 
    if (debugFlag) { // debugging... 
      double potN = potentialUnits*potential.get_val(locationN); 
      double potNEl = potentialUnitsElectr*electronPotential.get_val(locationN); 
      std::cerr << "---------------" << std::endl << " Decided on a step, start, vect6D ";
      for (int k=0; k !=6; ++k) std::cerr << ",  " <<  theVect6D[k];
      std::cerr << std::endl << " LocationS ";
      for (int k=0; k !=6; ++k) std::cerr << ",  " <<  locationS[k];
       std::cerr << std::endl << " .... Potential Here, Proton Bunch ref = "
	 << potS << " electron " << potSEl << " there " 
	 << potN << " electron " << potNEl <<  std::endl;
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

int RKIntegrator::propagateVElPy(PyObject *vect6DPy, 
                                 Real_scalar_field &phiBeam,
				 Real_scalar_field &phiElectron, 
                                 double tStart, double deltaT, PyObject *timePy) {
   double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
   double *tFinal = reinterpret_cast<double*>(PyArray_DATA(timePy));		 
   this->propagateV(loc, phiBeam, phiElectron, tStart, deltaT, tFinal);		 
}

int RKIntegrator::propThroughBunch(double *vect6D,  
                                   const Real_scalar_field &potential, 
                                   const Real_scalar_field &electronPotential, 
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
    nStep += this->propagateV(vect6D,  potential, 
                              electronPotential, tStart, 5.0e-10, &tFTmp);
    if (!inBeamPipe) break;
    tStart = this->getTime();
    kStep++;
  }
  (*tF) = tFTmp;
  return nStep;		 
}

int RKIntegrator::propThroughBunchWElPy(PyObject *vect6DPy, 
                                       Real_scalar_field &phiBeam,
				       Real_scalar_field &phiElectron,
                                       double tOffset, PyObject *tFPy) {
  		 
  double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
//  double rIn = std::sqrt(loc[0]*loc[0] + loc[2]*loc[2]); 
//  std::cerr << " RKIntegrator::propThroughBunch, rIn, C++ " << rIn << std::endl;   
   double *tFinal = reinterpret_cast<double*>(PyArray_DATA(tFPy));
   int nn = this->propThroughBunch(loc, phiBeam, phiElectron, tOffset, tFinal);
   return nn;		 
}

int RKIntegrator::propBetweenBunches(double *vect6D, const Real_scalar_field &phiElectron,
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
    double locationS[3], locationN[3];
    // Compute the field from the electron potential
    for (int k=0; k !=3; k++) { 
          locationS[k] = theVect6D[2*k];
          locationN[k] = theVect6D[2*k] + theVect6D[2*k+1]*speedOfLight*dt; 
    }  
    double potHere = potentialUnitsElectr*phiElectron.get_val(locationS);
    double potThere = potentialUnitsElectr*phiElectron.get_val(locationN);
    updateBField();
    for (int k=0; k !=3; k++) fieldHere[k] = 0.;
    for (int k=0; k !=3; ++k) fieldHere[k+3] = BFieldFromModel[k];
    if ((std::abs(potHere) > 0.01) || (std::abs(potThere) > 0.01)) {
      for (int k=0; k !=3; k++) {
        fieldHere[k] += -1.*eFieldUnitsElectr*phiElectron.get_deriv(locationS, k);
      }
    }
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
//      std::cerr << " prop between bunches Greater than 10000 integration steps.. " << std::endl;
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

int RKIntegrator::propBetweenBunchesElPy(PyObject *vect6DPy, 
                                         Real_scalar_field &electronPotential,
                                         double tOffset, double maxTime, 
					 PyObject *tFPy) {
		  
   double *loc = reinterpret_cast<double*>(PyArray_DATA(vect6DPy));
   double *tFinal = reinterpret_cast<double*>(PyArray_DATA(tFPy));		 
   return this->propBetweenBunches(loc, electronPotential, 
                                   tOffset, maxTime, tFinal);		  
}
void RKIntegrator::setUnitsElectrons(double totalCharge, double units0) {

// Total charge x 1.0/(4.0*pi*eps0) (in SI units) 
  potentialUnitsElectr = totalCharge*8.98755e9*4.0*M_PI; // Check with ../StudyPotential 
  eFieldUnitsElectr = potentialUnitsElectr/units0;
//  std::cerr << " Potential Units " << potentialUnits << " And quit " <<  std::endl;
//  exit(2); 
}
