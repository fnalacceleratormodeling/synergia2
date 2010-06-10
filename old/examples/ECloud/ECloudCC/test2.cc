#include <iostream>
#include <fstream>
#include <cmath>
#include "scalar_field.h"
#include "biGaussFunc.h"
//
// Testing the creation of a Real Scalar field, and numerical integration from 
// GSL.. By creating at two surimposed trigaussian.. 
// Since third dimension much bigger, set it to just a Gaussian along Z 


int main(int argc, char *argv[])  {

  double sig1[3]; double sig2[3]; double offset[3]; 
  sig1[0] = 1.5e-3;
  sig1[1] = 0.75e-3;
  sig1[2] = 1.5;
  sig2[0] = 1.5e-3;
  sig2[1] = 1.0e-3;
  sig2[2] = 0.75;
  //
  // Charge density, volumetric.. Irrelevant...  
  //
   
//  double dens1 = 8.9875e9*1.602e-19 *
//                 0.75e10 * 3.0/(4.0*M_PI*sig1[0]*sig1[1]*sig1[2]); 
//  double dens2 = 8.9875e9*1.602e-19 *
//                 1.5e10 * 3.0/(4.0*M_PI*sig2[0]*sig2[1]*sig2[2]); 
//  Using line density.. 
  double dens1 = 8.9875e9*1.602e-19 *
                 0.75e10 /(2.0*sig1[2]); 
  double dens2 = 8.9875e9*1.602e-19 *
                 1.5e10 /(2.0*sig2[2]); 
  //
  double gridSize[3];
  bool doFine=false;
  if (!doFine) { 
    for (int k=0; k != 3; k++) 
                   gridSize[k] = 10. * std::max(sig1[k], sig2[k]);
// Fine grid, to establish dependance at small radius .. 
  } else { 		   
    for (int k=0; k != 3; k++) 
                   gridSize[k] = 0.5 * std::max(sig1[k], sig2[k]);
  }		   
  double offsetPhi[3]; for (int k=0; k != 3; k++) offsetPhi[k] = 0.;		    
  for (int k=0; k != 3; k++) offset[k] = offsetPhi[k] - gridSize[k]/2.;
  int numPts[3]; numPts[0] = 400; numPts[1] = 400; numPts[2] = 256;
  std::cerr << " Test2, Grid Size " << gridSize[0] << ", " << 
                 gridSize[1] << " " << gridSize[2] << std::endl;
  Real_scalar_field myPhi(numPts, gridSize, offsetPhi);
  
  // Now compute the potential.
  int indices[3];
  biGaussFunc aBiG1(sig1, dens1);
  std::cerr << " First BiGaussian Integrals defined... " << std::endl;
  biGaussFunc aBiG2(sig2, dens2);
  std::cerr << " Second BiGaussian Integrals defined... " << std::endl;
  
  for (int iz=0; iz!=numPts[2]; iz++) {
    indices[2] = iz;
    double zz = offset[2] + iz*gridSize[2]/numPts[2];
    // Assume Gaussian.. Not likely to be right, but... 
    double factz1 = std::exp(-1.0*(zz*zz)/(2.0*sig1[2]*sig1[2]));
    double factz2 = std::exp(-1.0*(zz*zz)/(2.0*sig2[2]*sig2[2]));
    std::cerr << " At z slice " << iz << std::endl;
    for (int iy=0; iy!=numPts[1]; iy++) {
      indices[1] = iy;
      double yy = offset[1] + iy*gridSize[1]/numPts[1];
      for (int ix=0; ix!=numPts[0]; ix++) {
        indices[0] = ix;
        double xx = offset[0] + ix*gridSize[0]/numPts[0];
	double v1xy = aBiG1.val(xx, yy);
	double v2xy = aBiG2.val(xx, yy);
	double fieldVal = factz1*v1xy + factz2*v2xy;
//	std::cerr << " xx " << xx << " yy " << yy << " fieldVal " << fieldVal << std::endl;
	myPhi.get_points().set(indices, fieldVal);
      } // on ix 
    } // on iy
  } // On iZ 
  
  // Write it to a file... 
  std::string aFName("BiGaussianFieldV2.txt");
  myPhi.write_to_file(aFName);
  
  std::string fNameOut("AScalarFieldProfile");
  if (doFine) fNameOut += std::string("Fine");
  else fNameOut += std::string("Coarse");
  std::string fNameOut1(fNameOut); fNameOut1 += std::string(".txt");
  std::ofstream fOutC(fNameOut1.c_str());
  fOutC << " i x y z value " << std::endl;
  double stepSize[3]; for (int k=0; k!= 3; k++) stepSize[k] = 0.9*gridSize[k]/499;
  for (int i=0; i!=500; i++) {
    double location[3]; 
    for (int k=0; k!= 3; k++) location[k] = -0.45*gridSize[k] + stepSize[k]*i;
    fOutC << " " << i << " " << location[0] << " " << location[1] 
	               << " " << location[2] << " " << myPhi.get_val(location)
		       << std::endl;			      
  }
  fOutC.close();
  std::string fNameOut2(fNameOut); fNameOut2 += std::string("Radial.txt");
  std::ofstream fOutD(fNameOut2.c_str());
  fOutD << " i x y z value Ex Ey Ez " << std::endl;
  for (int i=0; i!=500; i++) {
    double location[3]; 
    for (int k=0; k!= 3; k++) location[k] = -0.45*gridSize[k] + stepSize[k]*i;
    fOutD << " " << i << " " << location[0] << " " << location[1] 
	              << " " << location[2] << " " << myPhi.get_val(location);
		       
    for (int k=0; k !=3; k++) fOutD << " " << myPhi.get_deriv(location, k);	
    fOutD  << std::endl;			      
  }
  fOutD.close();
 
}
