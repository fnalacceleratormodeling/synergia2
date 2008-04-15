#ifndef HAVE_CYLINDRICAL_H
#define HAVE_CYLINDRICAL_H

#include "macro_bunch_store.h"
#include "array_nd/array_2d.h"
#include "array_nd/array_3d.h"
#include "field_domain.h"

void get_cylindrical_coords(Macro_bunch_store &mbs, Array_2d<double> &coords);
void deposit_charge_cic_cylindrical(const Field_domain &fdomain, 
    Array_3d<double > &rho , Macro_bunch_store& mbs,
    const Array_2d<double> &coords);
void solve_tridiag_nonsym(const Array_1d<double> &diag,
    const Array_1d<double> &abovediag,
    const Array_1d<double> &belowdiag,
    const Array_1d<double> &rhs,
    Array_1d<double> &x);
void solve_cylindrical_finite_periodic(const Field_domain &fdomain,
    Array_3d<double> &rho, Array_3d<double> &phi);

#endif // HAVE_CYLINDRICAL_H
