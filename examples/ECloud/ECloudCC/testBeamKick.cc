#include <iostream>
#include <fstream>
#include <cmath>
#include "RKIntegrators.h"
#include "scalar_field.h"

// From Panagiotis Spentzouris, '2007
// Work in SI units, painfully slow and stupid, but ultimately easier to debug.. 
// Keep non-relativistic, though...  
// Using now a more formal class. 

int  main (void)
{

  
  Real_scalar_field myPhi; 
  // Obtained from test2 
  std::string aFName("BiGaussianFieldV2.txt");
  myPhi.read_from_file(aFName);
  std::string fNameOut1("CheckBiGaussinaFieldV2_1.txt");
  std::ofstream fOutC1(fNameOut1.c_str());
  fOutC1 << " i x y z value Ex Ey Ez " << std::endl;
  double stepSize[3]; stepSize[0] = 2.7e-5; stepSize[1] = 1.8e-5; stepSize[2] = 2.7e-2;
  double locStart[3]; locStart[0] = -6.75e-3; locStart[1] = -4.5e-3; locStart[2] = -6.75;
  double ETmp1[3];
  for (int i=0; i!=500; i++) {
    double location[3]; 
    for (int k=0; k!= 3; k++) location[k] = locStart[k] + stepSize[k]*i;
    for (int k=0; k !=3; k++) ETmp1[k] = myPhi.get_deriv(location, k);
    fOutC1 << " " << i << " " << location[0] << " " << location[1] 
	               << " " << location[2] << " " << myPhi.get_val(location);
    for (int k=0; k !=3; k++) fOutC1 << " " <<  ETmp1[0]; 
    fOutC1   << std::endl;			      
  }
  fOutC1.close();
//  std::cerr << " fOutC1 is closed .... And quit for now ...." << std::endl;
//  exit(2); 
  
  RKIntegrator myRK;
   
  double tStart = -10.0e-9;
  double tBigStep = 5.0e-10;  // by 1/2 of a ns, macro step...  
  // initial condition
  double betax = 0.; 
  double betay = 0.00015; 
  double betaz = 0.00005; 
  double y[6] = { 0.003, betax , 0.0, betay , 0.0, betaz};
  double yInit[6]; for (int k=0; k!=6; k++) yInit[k] = y[k];
  
  
  myRK.setToRelativstic(false);
  myRK.setPrecisionStep(1.0e-10);
  myRK.setDynamicRelativistic(true);
  myRK.setTrajectoryFileName("TrajBeamKickBiGaussV1j.txt");
  myRK.setMaximumXBeamPipe(5.5e-2);
  myRK.setMaximumYBeamPipe(2.5e-2);
  
  double tFinal = 0.;
  double explicitInterface=false;
  double bFieldOn = true;
  int kStep = 0;
  int allSteps = 0;
  if (explicitInterface) { 
    double tFTmp = tStart;
    while (tFTmp < 10.0e-9) {
      std::cerr << " Macro-step number kStep " << kStep << std::endl;
      allSteps += myRK.propagateV(y,  myPhi, tStart, tBigStep, &tFTmp);
      if (myRK.reachedBeamPipe()) {
        std::cerr << " Done bigtime step..Ok, reached beampipe, Clock at end " 
              << myRK.getTime() << std::endl;
        break;
      } else { 
        std::cerr << " Done bigtime step.... Keep going , Clock is now " 
              << myRK.getTime() << " tFinal tmp " << tFTmp << std::endl;
      }
      tStart = myRK.getTime();
      tFinal = tStart;
      kStep++;	      
    }
    myRK.closeTrajectoryFile();	      
  } else {
    double tOff = tStart;
    if (bFieldOn) {
      myRK.setBFieldStaticCmp(0.1, 1); // 1. kGauss
      // and let boost the electron to avoid boring trapping ? 
      y[3] = -0.004;
      y[1] = 0.0020;
     }
    allSteps = myRK.propThroughBunch(y,  myPhi, tOff, &tFinal); 
    myRK.closeTrajectoryFile();	      
    if (bFieldOn && (!myRK.reachedBeamPipe())) {
      std::cerr << " Starting propagation in with no bunch.. at t= " << tFinal << std::endl;
      myRK.reOpenTrajectoryFile("TrajBeamKickBiGaussV1gNF.txt");
      allSteps += myRK.propBetweenBunches(y, 0., 1.0e-8, &tFinal); // for an extra 10 ns. 
      myRK.closeTrajectoryFile();	      
    }
  }
  std::cerr << " Done Beam Kick test in " << kStep << " macro-steps " 
            << "( "<< allSteps << ") steps, current time " << tFinal << std::endl;
  
  return 0;
}

