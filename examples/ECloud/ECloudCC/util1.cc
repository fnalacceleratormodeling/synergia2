#include <iostream>
#include <cmath>

#include "util1.h"
#include "txionpack.h"

 gsl_rng * _saved_rng = 0;
 const gsl_rng_type * _saved_rng_type = 0;

void setGSLRandom() {

    _saved_rng = gsl_rng_alloc(gsl_rng_ranlxd2);
    _saved_rng_type =gsl_rng_ranlxd2;
    gsl_rng_set(_saved_rng, ((long int) 4934646));

}
// This is in fact bias Landau, reject the value less then 0.10
// No negative or ridiculously small values... 
// In addition, make sure we don't excess the maximum kinetic energy
// 
double getLandauEnergyDist(double width, double eMax) {

  double al=-1;
  while ((al < 0.10) || (al*width > eMax)) { 
   al= gsl_ran_landau(_saved_rng);
  }
  return al*width;
}
//
  
double getTechXEnergyDist(double eProton, double eMax) {

  double energy_in=eProton;
  double energy_outEl[1],  energy_outIon[1];
  double angle_outEl[1],  angle_outIon[1];
  energy_outEl[0] = eMax+1.0e10; 
  while (energy_outEl[0] > eMax) { 
    get_newparticles(1, 5, &energy_in, energy_outIon, energy_outEl, 
                              angle_outIon, angle_outEl);
  }			      
  return energy_outEl[0];
}

