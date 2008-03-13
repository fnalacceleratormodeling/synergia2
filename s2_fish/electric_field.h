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

#endif // HAVE_ELECTRIC_FIELD_H
