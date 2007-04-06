// -*- C++ -*-
/*******************************************
** solver.h
** Contains:
** 
*******************************************/

#ifndef HAVE_SOLVER_H
#define HAVE_SOLVER_H

#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>

double
fft_tester(int nx, int ny, int nz);

Real_scalar_field
solver_fft_open(Real_scalar_field &rho, bool periodic_z);

Real_scalar_field
solver_fd_multigrid_open(Real_scalar_field &rho);

Real_scalar_field
calculate_E_n(Real_scalar_field &phi, int n);

void
apply_E_n_kick(Real_scalar_field &E, int n_axis, double tau,
	       Macro_bunch_store &mbs);

#endif // HAVE_SOLVER_H
