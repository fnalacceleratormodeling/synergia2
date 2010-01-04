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
#include "fftw_helper.h"
#include <iostream>

Real_scalar_field
solver_fftw_open(Real_scalar_field &rho, Fftw_helper &fftw, bool periodic_z,
    bool use_guards=true);

#endif // HAVE_SOLVER_FFTW_H
