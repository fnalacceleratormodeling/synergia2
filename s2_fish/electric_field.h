// -*- C++ -*-
#ifndef HAVE_ELECTRIC_FIELD_H
#define HAVE_ELECTRIC_FIELD_H

#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>
#include "fftw_helper.h"

Real_scalar_field
calculate_E_n(Real_scalar_field &phi, int n,Fftw_helper &fftwh, bool z_periodic);

void
apply_E_n_kick(Real_scalar_field &E, int n_axis, double tau,
               Macro_bunch_store &mbs);

void
full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs,Fftw_helper &fftwh, bool z_periodic);

void
transverse_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs, Fftw_helper &fftwh, bool z_periodic);

void
full_kick_version(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs, Fftw_helper &fftwh, bool z_periodic);

void apply_Efield_kick(const std::vector<Real_scalar_field> &E, double tau,
               Macro_bunch_store &mbs);


//void testrho(Real_scalar_field &rho, Real_scalar_field &phi, Macro_bunch_store &mbs, bool z_periodic);


#endif // HAVE_ELECTRIC_FIELD_H
