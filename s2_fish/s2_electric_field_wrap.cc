#include "electric_field.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;

void
rw_kick_wrap(   double zsize,
                object &bin_partition,
                object &zdensity,
                object &xmom, 
                object &ymom,
                double tau, 
                Macro_bunch_store &mbs,
                double pipe_radius,
                double pipe_conduct, 
                double cutoff_small_z, 
                object &wake_coeff, 
                double quad_wake_sum, 
                bool quad_wake)
{
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
    rw_kick( zsize, bin_part_array,zdensity_array,
	    xmom_array,ymom_array,tau,mbs,pipe_radius,
            pipe_conduct, cutoff_small_z, wake_coeff_array,  quad_wake_sum, quad_wake);
}

BOOST_PYTHON_MODULE(s2_electric_field)
{
    def("calculate_E_n", calculate_E_n);
    def("apply_E_n_kick", apply_E_n_kick);
    def("full_kick", full_kick);
    def("transverse_kick", transverse_kick);
    def("rw_kick",rw_kick_wrap);
}

