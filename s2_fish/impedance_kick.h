// -*- C++ -*-
#ifndef HAVE_IMPEDANCE_KICK_H
#define HAVE_IMPEDANCE_KICK_H


#include "macro_bunch_store.h"
#include <iostream>


void
rw_kick(        double zsize,
                Array_1d<int> &bin_partition,
                Array_1d<double> &zdensity,
                Array_1d<double> &xmom, 
                Array_1d<double> &ymom,
                double tau, 
                Macro_bunch_store &mbs,
                double wake_factor,
                double cutoff_small_z, Array_1d<double> &wake_coeff,
                double quad_wake_sum,  bool bool_quad_wake);



void
get_wake_factors(int num_slices, int icut,
Array_1d<double> &zdensity, Array_1d<double> &xmom, Array_1d<double> &ymom,
Array_1d<double> &dipole_x, Array_1d<double> &dipole_y, Array_1d<double> &quad, Array_1d<double> &l_monopole);


#endif // HAVE_ELECTRIC_FIELD_H
