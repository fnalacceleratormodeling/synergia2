// -*- C++ -*-
#ifndef HAVE_COMMUNICATE_H
#define HAVE_COMMUNICATE_H

#include "scalar_field.h"
#include "fftw_helper.h"

#include <iostream>

void
gather_rho(Real_scalar_field &rho, int upper_limit);

void
fill_guards(Real_scalar_field &phi, Fftw_helper &fftwh);

void
broadcast_E(Real_scalar_field &E, int i_lower, int i_upper);

void
allgather_phi(Real_scalar_field &local_phi, Real_scalar_field &full_phi);

void
gather_global_rho(Real_scalar_field &local_rho,Real_scalar_field &global_rho);

#endif // HAVE_COMMUNICATE_H
