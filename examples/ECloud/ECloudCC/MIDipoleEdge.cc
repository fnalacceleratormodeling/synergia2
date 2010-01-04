#include "MIDipoleEdge.h"
#include <iostream>

MIDipQuadEdge::MIDipQuadEdge():miDipE(), MIEnergy(8.0) {

}

double MIDipQuadEdge::getFieldY(const double *location) {
  double zLoc = location[2] - zOffDipEdge;
  if (zLoc < 0.) {
    zLoc *= -1.;
    double byMax=BDipoleOverMIEnergy*MIEnergy;
    return (byMax*miDipE.getFieldYFact(zLoc));
  } else {
    double  yFact=miDipE.getFieldYFact(zLoc); // Assume the same z-function.. 
    double  stMax=BQuadStrOverMIEnergy*MIEnergy;
    return (stMax*location[1]*yFact);
  }
}

double MIDipQuadEdge::getFieldX(const double *location) {
  double zLoc = location[2] - zOffDipEdge;
  if (zLoc < 0.) {
    return 3e-4; // funky earthly magent field... 
  } else {
    double  yFact=miDipE.getFieldYFact(zLoc);
    double  stMax=BQuadStrOverMIEnergy*MIEnergy;
    return (stMax*location[0]*yFact);
  }
}

double MIDipQuadEdge::getFieldYFact(double z) {
  return miDipE.getFieldYFact(z);
}

MIDipoleEdge::MIDipoleEdge(){

  if (initDone) return;

  for (int k=0; k != numPts; k++) zs[k]*= .0254; // in metric..
  zMax *= .0254;
  double byAv=0.;
  for (int k=numPts-10; k != numPts; k++) byAv+=vals[k]; 
  byAv/=10.;
  for (int k=0; k != numPts; k++) vals[k]/= byAv; // in metric..
  initDone=true;
}

double MIDipoleEdge::getFieldYFact(double z) {

  if (z < 0.) return vals[0];
  if (z > zMax) return 1.0;
  int iz = (int) (z*numPts/zMax);
  if ((iz == 0) && (z < zs[1])) return (vals[0] + z*(vals[1]-vals[0])/zs[1]);
  if (iz == numPts-1) return 1.;
  if ((z > zs[iz]) && (z < zs[iz+1])) {
    return (vals[iz] + (z - zs[iz])*(vals[iz+1]-vals[iz])/(zs[iz+1]- zs[iz]));
  } else if (z < zs[iz]) {
//     std::cerr << " z " << z << " iz " << iz << std::endl;
     while (iz != numPts) {  
       iz--;
//       std::cerr << " iz ... " << iz << std::endl; 
       if ((z > zs[iz]) && (z < zs[iz+1])) 
         return (vals[iz] + (z - zs[iz])*(vals[iz+1]-vals[iz])/(zs[iz+1]- zs[iz]));
     }
     exit(2);
     return 1.; // should not occur. 
  } else if (z > zs[iz+1]) {
     while (iz != numPts) {  
       iz++;
       if ((z > zs[iz]) && (z < zs[iz+1])) 
         return (vals[iz] + (z - zs[iz])*(vals[iz+1]-vals[iz])/(zs[iz+1]- zs[iz]));
     }
     return vals[0]; // should not occur. 
  } 
}

//Data obtained from Dave Harding, Aug. 11 
// http://wwwtsmtf.fnal.gov/~dgcw/usranalysis/Archive/IDA/IDA032-0/special/ptscan/pt_red_pnts.212872.IDA032-0
 
bool  MIDipoleEdge::initDone=false;

double MIDipoleEdge::zs[]={ 0.000, 0.999, 1.994, 2.996, 3.741, 3.995, 4.243, 4.491, 4.739, 4.994, 
                            5.242, 5.490, 5.741, 5.989, 6.240, 6.492, 6.740, 6.991, 7.239, 7.487,	      
                            7.738, 7.990, 8.238, 8.486, 8.734, 8.988, 9.236, 9.488, 9.736, 9.984, 
			    10.232, 10.486, 10.731, 10.983, 11.234, 11.482, 11.733, 11.981, 12.230, 
                            12.481, 12.729, 12.983, 13.231, 13.483, 13.731, 13.982, 14.230, 14.481,  
                            14.729, 14.981, 15.229, 15.480, 15.728, 15.979, 16.227, 16.479, 16.727,  
                            16.978, 17.226, 17.477, 17.725, 17.977, 18.225, 18.476, 18.724, 18.975,  
                            19.224, 19.475, 19.475, 20.473, 21.472, 22.468, 23.470, 24.468, 25.467,  
                            26.466, 27.464, 28.463, 29.462, 30.461, 31.459, 32.458, 33.460, 34.455,  
                            35.457, 36.453, 37.451, 38.450, 39.449, 40.448, 41.446, 42.445, 43.444,  
			    44.445, 45.441, 46.440, 47.442, 48.440, 49.436, 50.435}; 

double MIDipoleEdge::vals[]={ 4.510697e-03, 5.688781e-03, 7.333306e-03, 9.781510e-03, 1.250359e-02, 
                              1.371241e-02, 1.510661e-02, 1.672785e-02, 1.862073e-02, 2.089353e-02, 
                              2.359248e-02, 2.684958e-02, 3.083987e-02, 3.581011e-02, 4.178801e-02,  
                              4.944510e-02, 5.868383e-02, 6.990535e-02, 8.257031e-02, 9.669603e-02,
			      1.118407e-01, 1.286032e-01, 1.459683e-01, 1.652512e-01, 1.851679e-01,  
			      2.069378e-01, 2.295686e-01, 2.536129e-01, 2.796822e-01, 3.068619e-01,  
                              3.365603e-01, 3.692985e-01, 4.041964e-01, 4.440482e-01, 4.887591e-01,  
                              5.387361e-01, 5.948987e-01, 6.587146e-01, 7.294441e-01, 8.066362e-01,  
                              8.878096e-01, 9.698261e-01, 1.046790e+00, 1.117693e+00, 1.179130e+00,  
                              1.230525e+00, 1.269904e+00, 1.301965e+00, 1.326240e+00, 1.342206e+00,  
                              1.353515e+00, 1.357472e+00, 1.360266e+00, 1.360203e+00, 1.361569e+00, 
                              1.361257e+00, 1.360663e+00, 1.360701e+00, 1.361267e+00, 1.361925e+00,  
                              1.362193e+00, 1.361160e+00, 1.361620e+00, 1.362546e+00, 1.361685e+00,  
                              1.364045e+00, 1.364661e+00, 1.364159e+00, 1.364471e+00, 1.365957e+00,  
                              1.367149e+00, 1.366340e+00, 1.366202e+00, 1.366693e+00, 1.367503e+00,  
                              1.367615e+00, 1.366953e+00, 1.368984e+00, 1.368043e+00, 1.366883e+00,  
                              1.367951e+00, 1.367494e+00, 1.367478e+00, 1.365422e+00, 1.367790e+00,  
                              1.366896e+00, 1.366292e+00, 1.367264e+00, 1.367168e+00, 1.366955e+00,  
                              1.367629e+00, 1.367660e+00, 1.367318e+00, 1.367326e+00, 1.367716e+00,
			      1.366566e+00, 1.366729e+00, 1.366376e+00, 1.365406e+00, 1.365736e+00};

double MIDipoleEdge::zMax=zs[MIDipoleEdge::numPts-1];			      

double MIDipQuadEdge::BDipoleOverMIEnergy=0.011409; // .09127/8.0, following Lattice file 
double MIDipQuadEdge::BQuadStrOverMIEnergy=-.132061; // -1.05649/8.0, following Lattice file 
// The following quantity has been guess-estimated the following way:
// A. The zip track data ( only a small fraction is defined above) shows that the 1/2 field 
// region extends from ziptrack coordinate 0.306 m and 6.382 m. Thus, the effective length 
// of the 6 m. dipole magnet is 6.076 m.  Lattice file says it is 6.09 m. Consistent.. 
// Thus, the magent center in ziptrack is at z=3.344 meter. 
// B. The typical distance between magents is ~0.7 m, edge to edge. Consistent with picture
//    Dave Capista showed me.  The legnth in the lattice file is 6.09 m. The center of the 
// of the magnet is therefore at 3.395 m. from our "in between magnet" point. 
// Thus, one has to shift only by 5.1 cm  
double MIDipQuadEdge::zOffDipEdge=0.051; // Don't know what the real value is yet! 
