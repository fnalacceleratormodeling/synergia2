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

Real_scalar_field
fft_tester(Real_scalar_field rho);

#endif // HAVE_SOLVER_H
