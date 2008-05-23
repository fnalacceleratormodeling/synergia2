#ifndef HAVE_ECLOUDCC_UTIL1_H
#define HAVE_ECLOUDCC_UTIL1_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>
//
// A set of utlity to create electrons, from the Gas and from the 
// walls. 
// April 2008


void setGSLRandom(); // Initialize my GSL Random number. 
double getLandauEnergyDist(double width, double eMax);
double getTechXEnergyDist(double eProton, double eMax);


#endif
