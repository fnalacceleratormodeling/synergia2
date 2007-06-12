#include <stdio.h>
#include "macro_bunch_store.h"

#include "BasErs_field.h"

using namespace boost::python;

int
apply_BasErs_kick(Macro_bunch_store &mbs, double sigmaX, double sigmaY, double sigmaZ, double gamma, double tau, double PartPersigmaZ)

{

// Line density will be calculated from sigma_z and 
// PartPersigmaZ which is the number of 
// charged particles per sigmaZ (signed for plus or minus 
// charges)

// sigmaX,Y,Z [m], at the accelerator frame

// Particles [internal units], convert here

  const double epsilon_0 = 8.854187817e-12;
  const double e_charge = 1.6021773e-19;

  //  double *get_funky_data = reinterpret_cast<double *>
  //  (reinterpret_cast<PyArrayObject*>(particles.ptr())->data);
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

  for (int ipart = 0; ipart < mbs.local_num; ipart++){
    for (int n_axis = 0; n_axis < 2; n_axis++){
      int index = 2 * n_axis + 1; // for axis n_axis = (0,1,2) corresponding to x,y,z,
      // in particle store indexing, px,py,pz = (1,3,5)
    
      std::vector<double> Efield(3);

      double LengthScale = mbs.units(0);

      xMeters = mbs.local_particles(0,ipart)/LengthScale;
      yMeters = mbs.local_particles(2,ipart)/LengthScale;
      Efield = myfield->NormalizedEField(xMeters, yMeters);

      std::cout << " Ex = " << Efield[0] << " Y = "<< Efield[1] << " Z = " << Efield[2] << std::endl; 

      //here we need to kick the friging y[coord]
      // we need the Efield back in the accelerator frame E->E*gamma
      // and a normalization of perv(2*scale)
      // Here is the explanation:
      // Leo normalizes E with  l/(2 pi epsilon_o), where l is the charge density C/m in the rest frame
      // Rob normalizes the scalar potential with ll/(4 pi epsilon_o), where ll is charge density 
      // at the lab frame

      //y[1] = y[1] + Efield[0]*Norm_field*tau*gamma;
      //y[3] = y[3] + Efield[1]*Norm_field*tau*gamma;

      double kick = Efield[n_axis]*Norm_field*tau*gamma;
      // update the data structure
      mbs.local_particles(index, n) += kick;
    }

  }

  delete myfield;
  
  return 0;
}

BOOST_PYTHON_MODULE(GaussSC)
{
  def("apply_BasErs_kick",&apply_BasErs_kick);
}
