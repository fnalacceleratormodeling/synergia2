#include <iostream>
#include <fstream>
#include <cmath>
#include "VorpalYeeField.h"
//
// Testing the Vorpal Yee Field interface. Just read a field..
// file name hardtype for now...and minimal test just make a simple table of the field 
// vx y of Ex, Ey, Ez  
//

int main(int argc, char *argv[])  {

  std::string fName("/local/lebrun/Synergia/VorpalStudies/myTest2/s1Es1ElM_YeeStaticElecFld_2.h5");
  
  VorpalYeeField aEMField(fName.c_str(), true, true);
  std::cerr << " Created VorpalYeeField ... " << std::endl;
  std::vector<double> physSize = aEMField.getPhysSize();
  std::vector<double> physOffset = aEMField.getPhysOffset();
  std::cerr << " Physical Offset, Y  " << physOffset[1] << " size " << physSize[1] 
            << " dim physOffset " <<physOffset.size() << std::endl;
  std::vector<double> loc(3);
  loc[0] = physOffset[0];
  loc[2] = physOffset[2];
  int nPts = 100;
  std::ofstream fOut("test3_YeeStaticElecFld_2.txt");
  fOut << " y Ex Ey Ez " << std::endl; 
  for (int k=0; k != nPts; k++) {
    loc[1] = physOffset[1] - physSize[1]/2. + k * (physSize[1]/nPts);
    fOut << " " << loc[1] << " " << aEMField.getEFieldX(loc) 
         << " " << aEMField.getEFieldY(loc) << " " << aEMField.getEFieldZ(loc) << std::endl; 
  }
  fOut.close();  
}
