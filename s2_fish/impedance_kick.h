// -*- C++ -*-
#ifndef HAVE_IMPEDANCE_KICK_H
#define HAVE_IMPEDANCE_KICK_H


#include "macro_bunch_store.h"
#include "array_nd/array_3d.h"
#include <iostream>


void
rw_kick(        Array_1d<double> &dparameters,
                Array_1d<int> &bin_partition,
                Array_1d<double> &zdensity,
                Array_1d<double> &xmom, 
                Array_1d<double> &ymom,
                Macro_bunch_store &mbs,
                 Array_1d<double> &wake_coeff,
                bool bool_quad_wake,
                int bunch_num,
                Array_3d<double> &stored_means,
                Array_2d<int>  &stored_buckets,
                Array_2d<double>  &stored_bunchnp
       );



void
get_wake_factors(int num_slices, int icut,
Array_1d<double> &zdensity, Array_1d<double> &xmom, Array_1d<double> &ymom,
Array_1d<double> &dipole_x, Array_1d<double> &dipole_y, Array_1d<double> &quad, Array_1d<double> &l_monopole,
 Array_3d<double> &stored_means, Array_2d<int>  &stored_buckets, Array_2d<double>  &stored_bunchnp,
 double line_length, double rbunch_sp,int bunch_i, double N_factor, double rescaling_f
                );


#endif // HAVE_ELECTRIC_FIELD_H
