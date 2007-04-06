#include <stdio.h>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <numarray/numarray.h>
#include <iostream>

#include "BasErs_field.h"

using namespace boost::python;

int
apply_BasErs_kick(numeric::array& particles, int numpart, double sigmaX, double sigmaY, double sigmaZ, double gamma, double tau, double PartPersigmaZ, double LengthScale)

{

// Line density will be calculated from sigma_z and 
// PartPersigmaZ which is the number of 
// charged particles per sigmaZ (signed for plus or minus 
// charges)

// sigmaX,Y,Z [m], at the accelerator frame

// Particles [internal units], convert here

  const double epsilon_0 = 8.854187817e-12;
  const double e_charge = 1.6021773e-19;

  double *get_funky_data = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(particles.ptr())->data);
  // for (i,j) element index = i*end_outer + j, with j 0 to end_outer

  double sigma[2];
  sigma[0]=sigmaX;
  sigma[1]=sigmaY;
  BasErs_field *myfield = new BasErs_field(sigma);

  // LineDensity at the beam frame
  double const LineDensity = PartPersigmaZ*e_charge/(sigmaZ*gamma);
  double const Norm_field = LineDensity/(M_TWOPI*epsilon_0);

  std::cout << " NormField = " << Norm_field << std::endl;

  double y[6];
  int N_PhS = 6; // Number of phase space coordinates

  for (int ipart = 0; ipart < numpart; ipart++){
    for (int j_PhS = 0; j_PhS < N_PhS; j_PhS++){
      y[j_PhS]=get_funky_data[ipart*N_PhS+j_PhS];
    }
    std::vector<double> Efield(3);

    xMeters = y[0]/LengthScale;
    yMeters = y[2]*LengthScale;
    Efield = myfield->NormalizedEField(xMeters, yMeters);

    std::cout << " Ex = " << Efield[0] << " Y = "<< Efield[1] << " Z = " << Efield[2] << std::endl; 

    //here we need to kick the friging y[coord]
    // we need the Efield back in the accelerator frame E->E*gamma
    // and a normalization of perv(2*scale)
    // Here is the explanation:
    // Leo normalizes E with  l/(2 pi epsilon_o), where l is the charge density C/m in the rest frame
    // Rob normalizes the scalar potential with ll/(4 pi epsilon_o), where ll is charge density 
    // at the lab frame

    y[1] = y[1] + Efield[0]*Norm_field*tau*gamma;
    y[3] = y[3] + Efield[1]*Norm_field*tau*gamma;

    // update the data structure
    for (int j_PhS = 0; j_PhS < N_PhS; j_PhS++){
      get_funky_data[ipart*N_PhS+j_PhS] =  y[j_PhS];
    }


  }

  delete myfield;
  
  return 0;
}

BOOST_PYTHON_MODULE(GaussSC)
{
  def("apply_BasErs_kick",&apply_BasErs_kick);
}
