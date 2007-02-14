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

Real_scalar_field
solver_fft_open(Real_scalar_field rho);

double
fft_tester(int nx, int ny, int nz);

#endif // HAVE_SOLVER_H
