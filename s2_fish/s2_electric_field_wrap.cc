#include "electric_field.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;

void
rw_kick_wrap(double zleft, double zsize,
                object &zdensity,
                object &xmom, 
                object &ymom,
                double tau, 
                Macro_bunch_store &mbs,
                double pipe_radiusx,
                double pipe_radiusy,
                double pipe_conduct,
                double zoffset)
{
   Array_1d<double> zdensity_array = 
        Array_nd_from_PyObject<double>(zdensity.ptr());
    Array_1d<double> xmom_array = 
        Array_nd_from_PyObject<double>(xmom.ptr());
    Array_1d<double> ymom_array = 
        Array_nd_from_PyObject<double>(ymom.ptr());
    rw_kick(zleft,zsize,zdensity_array,
	    xmom_array,ymom_array,tau,mbs,pipe_radiusx,
        pipe_radiusx,pipe_conduct,zoffset);
}

BOOST_PYTHON_MODULE(s2_electric_field)
{
    def("calculate_E_n", calculate_E_n);
    def("apply_E_n_kick", apply_E_n_kick);
    def("full_kick", full_kick);
    def("rw_kick",rw_kick_wrap);
}

