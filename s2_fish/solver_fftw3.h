// -*- C++ -*-
/*******************************************
** solver_fftw3.h
** Contains:
** 
*******************************************/

#ifndef HAVE_SOLVER_FFTW3_H
#define HAVE_SOLVER_FFTW3_H

#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>

namespace solver_fftw3 {
  Real_scalar_field
  solver_fftw3_open(Real_scalar_field &rho, bool periodic_z);
  void
  full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs);
  double
  fft_tester(int nx, int ny, int nz);
}

#endif // HAVE_SOLVER_FFTW3_H
