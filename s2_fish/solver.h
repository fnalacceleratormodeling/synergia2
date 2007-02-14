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
solver(Real_scalar_field rho);

double
fft_tester(int nx, int ny, int nz);

Real_scalar_field
fft_tester2(Real_scalar_field rho);

#endif // HAVE_SOLVER_H
