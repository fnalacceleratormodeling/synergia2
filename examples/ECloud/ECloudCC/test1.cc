#include <iostream>
#include <fstream>
#include <cmath>

#include "util1.h"
#include "txionpack.h"
int main(int argc, char *argv[])  {


  setGSLRandom();
  double meanL = 0.;
  double s2L = 0.;
  std::ofstream fOut("aLandauDistLow.txt");
  fOut << " i Landau TechX " << std::endl;
  int nSim = 100000;
  int nReal = 0;
  for (int i=0; i!=nSim; i++) {
     double al = getLandauEnergyDist(3., 500.);
     double eTx = getTechXEnergyDist(20.0e9, 500.);
     fOut << " " << i << " " << al << " " << eTx << std::endl;
     meanL += al;
     s2L += al*al;
     nReal++;
  }
  meanL /= nReal;
  double ss = std::sqrt((s2L - nReal*meanL*meanL)/(nReal-1.));
  std::cout << " mean " << meanL << " sigma " << ss << std::endl;
  fOut.close();
  // Try again, high energy tail...
  nSim = 1000000;
  int nGt2KevL = 0;
  for (int i=0; i!=nSim; i++) {
     double al = getLandauEnergyDist(5., 8.e7);
     if (al > 2.0e3) nGt2KevL++;
  }
  int nGt2KevT = 0;
  for (int i=0; i!=nSim; i++) {
     double eTx = getTechXEnergyDist(20.0e9, 8.7e7);
     if (eTx > 2.0e3) nGt2KevT++;
  }
  std::cout << " Number above 2 KeV, Landau  " << nGt2KevL 
       <<   " Tech-X " << nGt2KevT << std::endl;
       
  // Test angular/energy correlations with Tech-X... 
  nSim = 100000;
  std::ofstream fOutC("aTechXDist.txt");
  fOutC << " i e a " << std::endl;
  double eProton=2e10; // 20 GeV
  double energy_outEl[1],  energy_outIon[1];
  double angle_outEl[1],  angle_outIon[1];
  for (int i=0; i!=nSim; i++) {
    energy_outEl[0] = 1.0e10; 
    while (energy_outEl[0] > 8.7e7) { 
      get_newparticles(1, 5, &eProton, energy_outIon, energy_outEl, 
                              angle_outIon, angle_outEl);
    }
    fOutC << " " << i << " " << 
           energy_outEl[0] << " " << angle_outEl[0] << std::endl;			      
  }
  fOutC.close();
  
}
