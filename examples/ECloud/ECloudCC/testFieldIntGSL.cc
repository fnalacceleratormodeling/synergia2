#include <iostream>
#include <fstream>
#include <cmath>
#include "RKIntegrators.h"

// From Panagiotis Spentzouris, '2007
// Work in SI units, painfully slow and stupid, but ultimately easier to debug.. 
// Keep non-relativistic, though...  
// Using now a more formal class. 

int  main (void)
{

  // Debugging destructor ..

  RKIntegrator myRK;
   
//  field only... 

  double mu[6] = {0., 0., 0., 0., 0.0, 0.1}; // Bz = 1kG Gauss (0.1 T.) x sqrt(2) 

  double t = 0.0, t1 = 1.0e-7;  // 100 ns 
  double h = 1e-6; // 1 micron step size at the beginning.. 
  double betax = 0.; 
  double betay = 0.0015; 
  double betaz = 0.0005; 
  // Initial condition.. 
  double y[6] = { 0.003, betax , 0.0, betay , 0.0, betaz};
  double yInit[6]; for (int k=0; k!=6; k++) yInit[k] = y[k];
  
  
  myRK.setToRelativstic(false);
  myRK.setDynamicRelativistic(false);
  myRK.setTrajectoryFileName("GSLRKTestCstField_2b.txt");
  myRK.propagateF(y, mu, 1.0e-7);
  
  std::cerr << " Done first test.. " << std::endl;
 // Now run at 30 degree... 
  
  double mu2[6] = {0., 0., 0., 0.15, 0.866025*0.3, 0.0}; // Bz = 3kG Gauss (0.1 T.) 
  for (int k=0; k!=6; k++) y[k] = yInit[k];
  myRK.setTrajectoryFileName("GSLRKTestCst30degreeField_2b.txt");
  myRK.propagateF(y, mu2, 1.0e-7);
  
  std::cerr << " Done 2nd test.. " << std::endl;
  
  // Now accelerate and bend... B and E are parallel
  
  double mu3[6] = {0., 0., -1.0e4, 0.0, 0.0, 0.3}; // Bz = 3kG Gauss (0.3 T.), 10kV/m^2
  for (int k=0; k!=6; k++) y[k] = yInit[k];
  myRK.setTrajectoryFileName("GSLRKTestCstAccelField_2b.txt");
  myRK.propagateF(y, mu3, 2.0e-9);
  std::cerr << " Done 3rd test.. " << std::endl;
  
   // Now accelerate and bend... B and E are parallel
  
  double mu4[6] = {0., 0., -1.0e7, 0.0, 0.0, 0.3}; // Bz = 3kG Gauss (0.3 T.), 10MV/m^2
  for (int k=0; k!=6; k++) y[k] = yInit[k];
  myRK.setToRelativstic(true);
  myRK.setDynamicRelativistic(true);
  myRK.setTrajectoryFileName("GSLRKTestCstAccelRelField_10b.txt");
  myRK.resetClock();
  myRK.setPrecisionStep(1.0e-10);
  myRK.propagateF(y, mu4, 8.0e-10);
  std::cerr << " Done 4rd test..Clock at end " << myRK.getTime() << std::endl;
 
   // Now accelerate and bend... B vertical, E mostly are parallel
  
  double mu5[6] = {0., -1.0e7, 0.25e7, 0.0, 0.3, 0.}; // Bz = 3kG Gauss (0.3 T.), 10MV/m^2
  for (int k=0; k!=6; k++) y[k] = yInit[k];
  myRK.setToRelativstic(true);
  myRK.setDynamicRelativistic(true);
  myRK.setMaximumXBeamPipe(0.055);
  myRK.setMaximumYBeamPipe(0.025);
  myRK.setTrajectoryFileName("GSLRKTestCstAccelRelFieldBeamPipeCut_1c.txt");
  myRK.resetClock();
  myRK.setPrecisionStep(1.0e-10);
  myRK.propagateF(y, mu5, 8.0e-10);
  if (myRK.reachedBeamPipe())
    std::cerr << " Done 5rd test..Ok, reached beampipe, Clock at end " 
              << myRK.getTime() << std::endl;
    else 
     std::cerr << " Done 5rd test..?? No beam Pipe??? , Clock at end " 
              << myRK.getTime() << std::endl;
//  std::cerr << " Debugging Destructor.. " << std::endl;
//  exit(2);
  return 0;
}

