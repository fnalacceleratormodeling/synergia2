// -*- C++ -*-
/*******************************************
** solver_fftw2.h
** Contains:
** 
*******************************************/

#ifndef HAVE_SOLVER_FFTW2_H
#define HAVE_SOLVER_FFTW2_H

#include "scalar_field.h"
#include <iostream>

Real_scalar_field
solver_fftw2_open(Real_scalar_field &rho, bool periodic_z);

#endif // HAVE_SOLVER_FFTW2_H
