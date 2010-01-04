#ifndef HAVE_ECLOUD_RKINTEGRATORS_H
#define HAVE_ECLOUD_RKINTEGRATORS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/python.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "scalar_field.h"
#include "MIDipoleEdge.h"
#include "VorpalYeeField.h"

// Propagate electrons, using GSL integrator, for a short step. 

class RKIntegrator {
  public: 
    RKIntegrator(bool isRelativistic=false);
    ~RKIntegrator();
    
    // Returns the number of steps
    // vect6D is x, betax, y, betay, z, betaz.  Units are SI, meters. 
    // Bad old fashion C-code: input and output are the same vector. Saving a copy.. 
    // Fields is Ex, Ey, Ez, Bx, By, Bz, in the reference frame of the electron. 
    // This field is assumed perfectly static. 
  
    // Throughout these function, vect6D is both input and output. 
    // (initial and final position in phase space) 
    // Returns the number of steps 
    // (Basically, the number of calls to GSL numerical integration 
    
    int propagateF(double *vect6D, const double *fields, double deltaT);
    // Same with Python interface. 
    int propagateFPy(PyObject *vect6DPy, PyObject *fieldsPy, double deltaT);
  
    // Propagate from the potential of the bunch at the point given where the electron is...
    // from an arbitrary time, which can be negative, for an approximate time deltaT. 
    // depending on accuracy requirement, tFinal (returned) is approximately tStart + tFinal. 
    // That is, if the electron did not bumped in the beam pipe before, at which point 
    // This method is based on propagateF 
  
    int propagateV(double *vect6D,  const Real_scalar_field &potential, 
                 double tStart, double deltaT, double *tFinal);
    int propagateVPy(PyObject *vect6DPy, Real_scalar_field &fieldsPy, 
                 double tStart, double deltaT, PyObject *timePy);
		 
    // Same as above, but adding the field from the electron 
    // As for the proton, the electron field is assumed 
  
    int propagateV(double *vect6D,  const Real_scalar_field &potential, 
                   const Real_scalar_field &electronPotential,
                 double tStart, double deltaT, double *tFinal);
    int propagateVElPy(PyObject *vect6DPy, Real_scalar_field &fieldsPy, 
                  Real_scalar_field &fieldsElectronPy, 
                 double tStart, double deltaT, PyObject *timePy);
    //
    // Propagate through the Proton bunch defined by a given potential
    // Start at tOffset, return tFinal.. We will step through... 
    // we will stop when we are outside the region where the protential 
    // is non-zero. 
    // This method is based on propagateV
    //
    int propThroughBunch(double *vect6D,  const Real_scalar_field &potential, 
                 double tOffset, double *tFinal);
    int propThroughBunchPy(PyObject *vect6DPy,  Real_scalar_field &potential, 
                 double tOffset, PyObject *tFinal);
    // And with Electrons 	   
    int propThroughBunch(double *vect6D,  const Real_scalar_field &potential, 
                  const Real_scalar_field &electronPotential, 
                 double tOffset, double *tFinal);
    int propThroughBunchWElPy(PyObject *vect6DPy,  Real_scalar_field &potential,
                 Real_scalar_field &electronPotential,
                 double tOffset, PyObject *tFinal);
   //		 
   // For comparing Vorpal vs this in terms of propagating accuracy.. 
   // 
   int propThroughVorpalField(double *vect6D, const VorpalYeeField &aField, 
                              double tStart, double deltaT, double *tFinal); 		 
   //
    // Propagate in between bunch, where there are no electric field, and 
    // the magnetic field is static, and uniform. bField is of type double[3] 
    //
    int propBetweenBunches(double *vect6D,   
                         double tOffset, double maxTime, double *tFinal);
    int propBetweenBunchesPy(PyObject *vect6DPy,  
                           double tOffset, double maxTime, PyObject *tFinal);
			   
   // And with electrons 		   
    int propBetweenBunches(double *vect6D, const Real_scalar_field &electronPotential,  
                         double tOffset, double maxTime, double *tFinal);
    int propBetweenBunchesElPy(PyObject *vect6DPy, Real_scalar_field &electronPotential, 
                           double tOffset, double maxTime, PyObject *tFinal);
			   
    enum BFieldModel {UNIFORM, DIPOLEEDGE, QUADRUPOLE, QUADRUPOLEEDGE, DIPQUADEDGES};			   
		   
    void closeTrajectoryFile(); 
    void reOpenTrajectoryFile(const char *fName); 
    void reOpenTrajectoryFilePy(PyObject *fName); 
    void setTrajectoryFileNamePy(PyObject *fName);	 
    PyObject *getTrajectoryFileNamePy();
    // The Synergia potential have natural units. 
    // These units are not available in the scalar_field class, 
    // so the have to be uploaded here..
    // Bunch charge here is expressed in number of protons..
    // units0 is the related to "frequency scale" in the Impact package.
    // The charge of the bunch must be expressed in Coulomb.
    // Units0 is a scale factor, related to Impact frequency_scale paramter.  	 
    void setUnits(double totalChargeProtonBunch, double units0); // 
    void setUnitsElectrons(double totalChargeElectron, double units0); // 
    void setFieldModel(BFieldModel aModel, double strength);
    // I don't know how to handle C++ enumeration in Python... Sorry... 
    void setFieldModelPy(int aModel, double strength);
 
  private:
    
    static const double speedOfLight; // m/sec
    static const double speedOfLightSq; // (m/sec)^2
    static const double electronCharge; // Coulomb
    static const double electronMass;
    static const double electronMassEV;
    static const double ChargeOverMass;
    static int signChange; // Flip the sign of electron default is 1 
    static const gsl_odeiv_step_type *myOdeiv_step_type;

    mutable bool isRel;  // Relativistic or not.  very dynamic.
    double relCutOff; 
    bool trajToFile;   // Dump the trajectory to a file 
    std::string trajFName;
    mutable double precisionStep; // set by the GSL algorithms. 
                                  // We keep here a value set by the user. 
    mutable double stepRatio; // Avoid jumping too fast in numrical integration, adjusting 
                              // the macro-step...Obsolete, replace by Laromor radius criteria. 
    mutable double minLarmorRad; // below this, we will not keep track of the Larmor phase 
                                 // anylonger ( it will be chosen randomly). 			      			  
    mutable bool inBeamPipe;
    mutable bool theDebugFlag; // Lots of print out...
    mutable bool errorInStep; // Bad nes, field boundary problems
    double gamProtonBunch;
    double gamProtonBunchSq;
    BFieldModel staticFieldModel;
    double BFieldStatic[3];
    double BFieldFromModel[3];
    double BFieldStrength; // Model dependant, gradient for quadrupole... 
    double BFieldMaxNonUnif; // to compute the maximum value of Larmor radius.
    double BFieldMinNonUnif; // and its minimum..
    double theTAbsolute;  // for integrating on dTau 
    double theTAbsoluteB;  // Boosted time, for above... 
    double theTStart; // for moving bunch...
    double theDeltaTGoal; // Maximum delta for macro step, using to terminate recursion 
                         // properly...
    double theTimePrecisionIntegralLimit; 			  
    mutable double theCPUTimeStart;
    mutable double theCPUTimeEnd;
    
    mutable double theCurrentPotential; // The potential in the fine RK step. Used for debugging... 
    double theVect6D[6]; // the vector at current time... 
                      // Changed in propagateStepField
    bool relDynamic; // flag to dynamically change the state Relativistic/non-relativistic.
    // Geometry constraint..
    double maxXBeamPipe;
    double maxYBeamPipe;
    double minMaxRPipe;
    int nRecursive; // for debugging only... 

    // More parameters to control accuracy..
    // Specifiy the tolerance for field uniformity over a presumed time step. 
    // Used in propagateV   
    double relEFieldChange;
    // Unit convesion factors. Synergia inherits from Impact, 
    // which uses "natural units" that depend on grid spacing, 
    // or non-standard units of charge, vaccum permeability, and what not.. 
    double potentialUnits; // For the potential
    double eFieldUnits; // for its derivative.. done O.K., should not be 
                        // needed.. 
// Same for the electrons. 
    double potentialUnitsElectr; // For the potential
    double eFieldUnitsElectr;
 
    std::ofstream *thefOutPtr;
    // Pointer to GSL utilities.. 
    gsl_odeiv_step * sOdeivStepAlloc; 
    gsl_odeiv_control * cOdeivControl; 
    gsl_odeiv_evolve * eOdeivEvolve; 
    // Magnetic field knowledge... Make it a pointer... 
    // Kind of ugly, this class is supposed to be end-application independent.. 
    MIDipQuadEdge *MIMagnetData; // Owned 
        
    mutable std::vector<double> currentField;
    		
    // Returns the number of steps. 	
    int propagateStepField(const double *fields, double deltaT);
    int propagateFastNL(const double *fields, double deltaT, double rho);
    // Interface to GSL Ordinary Differential equation Integrators 
    static int funcNR (double t, const double y[], double f[], void *params);
    static int funcR (double t, const double y[], double f[], void *params);
    static int jacNR (double t, const double y[], double *dfdy, 
                double dfdt[], void *params);
    static int jacR (double t, const double y[], double *dfdy, 
                double dfdt[], void *params);
    void updateBField(); // Based on model.. No longer uniform, but still static. 
    void fieldBunchTransBeamToLab(const double *fieldIn, double *fieldOut);
    bool checkBeamPipeBoundary(const double *vect6D);
    double distToBeamPipe(const double *vect6D);
    double getLarmorRadius(double *vect6D);  
    void writeFOutTrajPreamble();
    double getcputime(); // for improving integration speed.

  public:

  inline void setToRelativstic(bool GoRelativistic) { isRel=GoRelativistic; }   
  inline bool isRelativistic() const {return isRel;}
  inline bool getRelativisticCutOff() const {return relCutOff;}
  inline bool setRelativisticCutOff(double c) {relCutOff=c;}
  inline void setPrecisionStep(double h) {
      precisionStep = h;
      gsl_odeiv_control_free (cOdeivControl);
      cOdeivControl = gsl_odeiv_control_y_new (precisionStep, 0.0); 
  }
  //
  // Tweaking the accuracy with which the upper integration limit is reached... 
  // 
  
  inline void setTheTimePrecisionIntegralLimit( double deltaTLim) {theTimePrecisionIntegralLimit = deltaTLim;}
  inline double getTheTimePrecisionIntegralLimit() const { return theTimePrecisionIntegralLimit;}
  
  inline double getPrecisionStep() const {return precisionStep; }
  inline void setTrajectoryFileName(const char *fName) 
     {this->reOpenTrajectoryFile(fName); }
  inline const char* getTrajectoryFileName() const {return trajFName.c_str();}

  inline void setGamProtonBunch(double gamma) 
     {gamProtonBunch=gamma; gamProtonBunchSq = gamma*gamma;} 
  inline double getGamProtonBunch() const { return gamProtonBunch;} 

  inline void setBFieldStaticCmp(double val, int axis)
        { BFieldStatic[axis] = val;
	  double bNorm = 0.;
	  for (int k=0; k!=3; k++) bNorm += BFieldStatic[k]*BFieldStatic[k];
	  BFieldMaxNonUnif=std::sqrt(bNorm); 
	  BFieldMinNonUnif=BFieldMaxNonUnif; } // No protection,, Only valid for uniform Model
	   
  inline double getBFieldStaticCmp(int axis) const {return BFieldStatic[axis];}		

  inline double getStepRatio() const { return stepRatio;}
  inline void setStepRatio(double r) { stepRatio = r;}
  
  inline double getMinimumLarmorRadius() const { return minLarmorRad;}
  inline void setMinimumLarmorRadius(double r) { minLarmorRad = r;}
  
  // Allowing to change the equation set from non-relativistic to relativistic
  // dynamically, based on current velocity. 
  inline bool getDynamicRelativistic() const { return relDynamic;}
  inline void setDynamicRelativistic(bool goRelDynamic) { relDynamic = goRelDynamic;}
  // 
  // Dynamic specification of the step size 
  // by requesting a specific spatial uniformity. Over this time step, 
  // the maximum relative variation of any of the electric field componenent 
  // will be less than..
  inline double getRelEFieldChange() const {return relEFieldChange;}  
  inline void setRelEFieldChange(double rr) {relEFieldChange = rr;}
  
  
  inline double getMaximumXBeamPipe() const { return maxXBeamPipe;}
  inline double getMaximumYBeamPipe() const { return maxYBeamPipe;}
  
  inline void setMaximumXBeamPipe(double xm) { maxXBeamPipe=xm;}
  inline void setMaximumYBeamPipe(double ym) { maxYBeamPipe=ym;}
  
  inline bool reachedBeamPipe() const {return !inBeamPipe;}
  
  inline void resetClock() {theTAbsolute=0.;  theTAbsoluteB=0.;}
  inline double getTime() const {return theTAbsolute;}
  inline double getTimeBoosted() const {return theTAbsoluteB;}
  
  inline void setDebugOn() const {theDebugFlag=true;}
  inline void setDebugOff() const {theDebugFlag=false;}
  // For debugging from the static method...
  inline double getCurrentVect6D(int k) const {return theVect6D[k];}
  
  
  inline bool gotPropagationError() const {return errorInStep;}
  inline void setToPositron() {signChange = -1;}
  inline void setToElectron() {signChange = 1;}
  // scale factor for the electric field from the electrons
    
  inline BFieldModel getFieldModel() const {return staticFieldModel;}
  inline double getCPUTimeLastTrack() const {return (theCPUTimeEnd - theCPUTimeStart);}
};
#endif
