#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "VorpalElectrons.h"
#include "VorpalYeeField.h"
#include "RKIntegrators.h"
#include "util1.h"
//
// Testing the Vorpal Electron interface, and later, the comparison of the propagation 
// of electrons Vorpal/RKIntergator interface.
//

int main(int argc, char *argv[])  {

  setGSLRandom();
  std::string dirName("/local/lebrun/Synergia/VorpalStudies/myTest2/");
  std::string token("s1Es1ElMBp1Y2p25_1kV_"); // The Vorpal ".pre" name 

//  std::cout << " Looping electrons file Number ... " << std::endl;
// read the static voltage once..
  std::string fNameField(dirName); 
  fNameField += token + std::string("YeeStaticElecFld_1.h5");   
  VorpalYeeField aEMField(fNameField.c_str(), true, true);
//  std::vector<double> eFieldCst(3); eFieldCst[0] = -5000.; eFieldCst[1] = -5000.; eFieldCst[2] = 0.;
//  aEMField.addCstEField(eFieldCst);
  // Check EField at coordinate where disagreement occurs... 
  std::vector<double> loc(3); loc[0] = .0074; loc[1]=1.5e-3; loc[2]=1.0e-4;
  std::cout << " Electric Field "; 
  for (int k=0; k !=3; k++) std::cout << " " << aEMField.getEField(loc,k); 
  std::cout << std::endl;
  std::vector<double> physSize = aEMField.getPhysSize();
  std::cout << " Physical size of the field ";
  for (int k=0; k !=3; k++) std::cout << " " << physSize[k]; 
  std::cout << std::endl;
 
  // Integrate the dum way along y=0., y=0. 
  if(token.find("B0Y0") != std::string::npos) {
    double vXTot = 0.;
    double xx = 5.e-9;
    double dx = 1.0e-6;
    std::vector<double> loc(3); loc[0] = xx; loc[1]=1.5e-9; loc[2]=1.0e-9;
    while (xx < .0130043) {
      loc[0] = xx;
      vXTot += dx*aEMField.getEFieldX(loc);
      xx += dx;
    }
    std::cout << " Integrated voltage up to .0130043 " << vXTot << std::endl;
    vXTot = 0.;
    xx = 5.e-9;
    dx = 1.0e-5;
    while (xx < .01427910) {
      loc[0] = xx;
      vXTot += dx*aEMField.getEFieldX(loc);
      xx += dx;
    }
    std::cout << " Integrated voltage up to .0147910 " << vXTot << std::endl;
  }
  
  double vect6D[6];
  double t0=0.;
  
  
  RKIntegrator myRK;
  myRK.setDebugOff();
  myRK.setToRelativstic(true);
  myRK.setRelativisticCutOff(1.0e-14);
  myRK.setRelEFieldChange(5.0e-3);
//  double aRKFieldBgrd[6] = {0., 0., 0., .01, 0.0, 0.0}; // Bz = 100 Gauss (0.01 T.) 
  double aRKFieldBgrd[6] = {0., 0., 0., 0., 0.0, 0.0}; // Bz = 100 Gauss (0.01 T.) 
  if (token.find("Bp1Y") != std::string::npos) aRKFieldBgrd[4] = 0.1;
  for (int k=0; k!=3; k++) myRK.setBFieldStaticCmp(aRKFieldBgrd[k+3], k);
  myRK.setFieldModel(RKIntegrator::UNIFORM, std::sqrt(aRKFieldBgrd[0]*aRKFieldBgrd[0] 
                     + aRKFieldBgrd[1]*aRKFieldBgrd[1] + aRKFieldBgrd[2]*aRKFieldBgrd[2]));
  
  double tPrevious = 0.;
  double tPrevRK = 0.;
  double tRKAfterStep = 0.;
  bool modeTime = true; 
  if (modeTime) { 
    std::string fOutName("./test4_"); fOutName += token;
    fOutName += std::string("DynRel12CC.txt");
    std::ofstream fOut(fOutName.c_str());
    fOut << " t  XV   betaXV  YV   betaYV  ZV   betaZV   eKV  ";
    fOut << " tRK XS   betaXS   YS   betaYS  ZS   betaZS    eKS  ";
    fOut << std::endl;
    for (int kFile=1; kFile !=14; kFile++) { 
      std::ostringstream kfStr; kfStr << kFile;
      std::string fName (dirName); fName += token + std::string("electrons_") + kfStr.str() + std::string(".h5"); 
      VorpalElectrons ees(fName.c_str());
      if (!ees.isFileValid()) break;
      std::vector<aVorpalElectron>::const_iterator ie = ees.begin();
      fOut << " " << ees.getTime() << " " << ie->x << " " << ie->betax 
              << " " << ie->y << " " << ie->betay
	      << " " << ie->z << " " << ie->betaz << " " << ie->kineticEnergy;
     // Define initial condition. 
     if (kFile == 1) {
       const double *eDbl = &(ie->x);
       for (int k=0; k!=6; k++) vect6D[k] = eDbl[k]; 
       t0=ees.getTime(); 
       tPrevious = t0; tPrevRK = t0; tRKAfterStep=t0;
     } else {
     // Synergia (RKIntegrator)  propagate..
       double dt = ees.getTime() - tPrevious;
       double tTmp;
       myRK.propThroughVorpalField(vect6D, aEMField, 0., dt, &tTmp);
       tRKAfterStep += tTmp;
       std::cerr << " After propThroughVorpalField, macro step " << kFile << " betaX " 
               << vect6D[1] << std::endl;
//	if (vect6D[1] > .029) myRK.setDebugOn();       
     }
     fOut << " " << tRKAfterStep;
     for (int k=0; k != 6; k++) fOut << " " << vect6D[k];
     double betaSqS = 0.; 
     for (int k=0; k!=3; k++) {
       betaSqS += vect6D[2*k+1]*vect6D[2*k+1];
     }
     double gammaS = std::sqrt(1.0/(1.0-betaSqS));
     double eKS = 510998.918 * (gammaS - 1.0); // need energy in GeV to compute rho.. 
     fOut << " " << eKS << std::endl;
     tPrevious = ees.getTime();
     tPrevRK = tRKAfterStep; 
    }// on step or time 
    fOut.close();
  } else {
    std::string fOutName("./test4_Pnum_"); fOutName += token + std::string(".txt");
    std::ofstream fOut(fOutName.c_str());
    fOut << " tpNum  t XV   betaXV  YV   betaYV  ZV   betaZV   eKV  ";
    fOut << " tRK XS   betaXS   YS   betaYS  ZS   betaZS    eKS  ";
    fOut << std::endl;
    int kMaxFile=3;
    for (int iPart=0; iPart != 11; iPart++) {
      for (int kFile=1; kFile != (kMaxFile+1); kFile++) { 
        std::ostringstream kfStr; kfStr << kFile;
        std::string fName (dirName); fName += token + std::string("electrons_") + kfStr.str() + std::string(".h5"); 
        VorpalElectrons ees(fName.c_str());
        if (!ees.isFileValid()) break;
        std::vector<aVorpalElectron>::const_iterator ie = ees.begin();
        if (iPart > 0) for (int kpp=0; kpp!=iPart; kpp++) ie++;
        if (kMaxFile == kFile) fOut << " " << iPart << " " << ees.getTime() << " " << ie->x << " " << ie->betax 
              << " " << ie->y << " " << ie->betay
	      << " " << ie->z << " " << ie->betaz << " " << ie->kineticEnergy;
       // Define initial condition. 
       if (kFile == 1) {
         const double *eDbl = &(ie->x);
         for (int k=0; k!=6; k++) vect6D[k] = eDbl[k]; 
         t0=ees.getTime(); 
         tPrevious = t0; tPrevRK = t0; tRKAfterStep=t0;
       } else {
       // Synergia (RKIntegrator)  propagate..
         double dt = ees.getTime() - tPrevious;
         double tTmp;
         myRK.propThroughVorpalField(vect6D, aEMField, 0., dt, &tTmp);
         tRKAfterStep += tTmp;
         std::cerr << " After propThroughVorpalField, macro step " << kFile << " betaX " 
               << vect6D[1] << std::endl;
         tPrevious = ees.getTime();
         tPrevRK = tRKAfterStep; 
//	if (vect6D[1] > .029) myRK.setDebugOn();       
       }
	
     } // on time (i.e., Vorpal output files. 
     fOut << " " << tRKAfterStep;
     for (int k=0; k != 6; k++) fOut << " " << vect6D[k];
     double betaSqS = 0.; 
     for (int k=0; k!=3; k++) {
       betaSqS += vect6D[2*k+1]*vect6D[2*k+1];
     }
     double gammaS = std::sqrt(1.0/(1.0-betaSqS));
     double eKS = 510998.918 * (gammaS - 1.0); // need energy in GeV to compute rho.. 
     fOut << " " << eKS << std::endl;
    } // on particles within that file 
  } // on time or space option
}
