#include "impedance_kick.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;

void
rw_kick_wrap(  
                object &dparameters,
                object &bin_partition,
                object &zdensity,
                object &xmom, 
                object &ymom,
                Macro_bunch_store &mbs,
                object &wake_coeff, 
                bool bool_quad_wake,
                int bunch_i,
                object &stored_means,
                object &stored_buckets,
                object &stored_bunchnp
            )   
{

   Array_1d<double> dparameters_array=Array_nd_from_PyObject<double>(dparameters.ptr());
   Array_1d<double> zdensity_array = 
        Array_nd_from_PyObject<double>(zdensity.ptr());
    Array_1d<double> xmom_array = 
        Array_nd_from_PyObject<double>(xmom.ptr());
    Array_1d<double> ymom_array = 
        Array_nd_from_PyObject<double>(ymom.ptr());
    Array_1d<int> bin_part_array = 
        Array_nd_from_PyObject<int>(bin_partition.ptr());
    Array_1d<double> wake_coeff_array = 
        Array_nd_from_PyObject<double>(wake_coeff.ptr());
    Array_3d<double> stored_means_array=
            Array_nd_from_PyObject<double>(stored_means.ptr());
    Array_2d<int> stored_buckets_array = 
            Array_nd_from_PyObject<int>(stored_buckets.ptr()); 
    Array_2d<double> stored_bunchnp_array=
            Array_nd_from_PyObject<double>(stored_bunchnp.ptr());  
    
    rw_kick( dparameters_array, bin_part_array,zdensity_array,
	  xmom_array,ymom_array,mbs,
      wake_coeff_array,   bool_quad_wake, bunch_i,
      stored_means_array, stored_buckets_array, stored_bunchnp_array);
    
//     rw_kick( zsize, bin_part_array,zdensity_array,
//              xmom_array,ymom_array,tau,mbs,wake_factor,
//              cutoff_small_z, wake_coeff_array,  quad_wake_sum, bool_quad_wake, bunch_num,
//              stored_means_array, stored_buckets_array,bunch_spacing);
}

BOOST_PYTHON_MODULE(s2_impedance_kick)
{
    def("rw_kick",rw_kick_wrap);
}

