// -*- C++ -*-
/*******************************************
** solver_fftw.h
** Contains:
** 
*******************************************/

#ifndef HAVE_SOLVER_FFTW_H
#define HAVE_SOLVER_FFTW_H

#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>

Real_scalar_field
solver_fftw_open(Real_scalar_field &rho, bool periodic_z);

#endif // HAVE_SOLVER_FFTW_H
