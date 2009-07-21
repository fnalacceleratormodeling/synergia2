// -*- C++ -*-
/*******************************************
** deposit.h
** Contains:
**
*******************************************/

#ifndef HAVE_DEPOSIT_H
#define HAVE_DEPOSIT_H

#include "array_nd/array_1d.h"
#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
double
deposit_charge_cic(Real_scalar_field& sf, Macro_bunch_store& mbs,
                   bool z_periodic);

// Deposit charge in the bin to the left hand side of the particle
double
deposit_charge_ngp(Real_scalar_field& sf, Macro_bunch_store& mbs);

void
rho_to_rwvars(Real_scalar_field &rho, Array_1d<double> &zdensity,
            Array_1d<double> &xmom, Array_1d<double> &ymom);

void
calculate_rwvars(Macro_bunch_store& mbs,
                   Array_1d<double> &zdensity,
                   Array_1d<double> &xmom, Array_1d<double> &ymom,
                   double z_left, double z_length, Array_1d<int> &slice_partition);

#endif // HAVE_DEPOSIT_H
