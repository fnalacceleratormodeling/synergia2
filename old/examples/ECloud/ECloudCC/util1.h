#ifndef HAVE_ECLOUDCC_UTIL1_H
#define HAVE_ECLOUDCC_UTIL1_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
//
// A set of utlity to create electrons, from the Gas and from the 
// walls. 
// April 2008

extern  gsl_rng * _saved_rng;

void setGSLRandom(); // Initialize my GSL Random number. 
double getLandauEnergyDist(double width, double eMax);
double getTechXEnergyDist(double eProton, double eMax);
inline double getGSLRanFlatPhi() {
  return gsl_ran_flat(_saved_rng, 0., 2.0*M_PI);
}
#endif
