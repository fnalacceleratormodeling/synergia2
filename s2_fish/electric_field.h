// -*- C++ -*-
#ifndef HAVE_ELECTRIC_FIELD_H
#define HAVE_ELECTRIC_FIELD_H

#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>

Real_scalar_field
calculate_E_n(Real_scalar_field &phi, int n);

void
apply_E_n_kick(Real_scalar_field &E, int n_axis, double tau,
               Macro_bunch_store &mbs);

void
full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs);

void
full_kick_version(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs);

void
rw_kick(double zleft, double zsize,
                Array_1d<double> &zdensity,
                Array_1d<double> &xmom, 
                Array_1d<double> &ymom,
                double tau, 
                Macro_bunch_store &mbs,
                double pipe_radiusx,
                double pipe_radiusy,
                double pipe_conduct,
                double zoffset);

void apply_Efield_kick(const std::vector<Real_scalar_field> &E, double tau,
               Macro_bunch_store &mbs);

#endif // HAVE_ELECTRIC_FIELD_H
