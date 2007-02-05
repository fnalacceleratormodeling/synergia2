/*******************************************
** deposit.h
** Contains:
** 
*******************************************/

#ifndef HAVE_DEPOSIT_H
#define HAVE_DEPOSIT_H

#include "scalar_field.h"
#include "macro_bunch_store.h"
#include <iostream>

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
double
deposit_charge_cic(Real_scalar_field& sf, Macro_bunch_store& mbs);

// Deposit charge in the bin to the left hand side of the particle
double
deposit_charge_ngp(Real_scalar_field& sf, Macro_bunch_store& mbs);

#endif // HAVE_DEPOSIT_H
