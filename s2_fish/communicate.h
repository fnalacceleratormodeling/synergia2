// -*- C++ -*-
#ifndef HAVE_COMMUNICATE_H
#define HAVE_COMMUNICATE_H

#include <iostream>
// undefine symbols that conflict between iostream and mpich2
#if defined(SEEK_CUR)
#undef SEEK_CUR
#undef SEEK_SET
#undef SEEK_END
#endif /* defined(SEEK_CUR) */

#include "fftw_helper.h"
#include "scalar_field.h"

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
