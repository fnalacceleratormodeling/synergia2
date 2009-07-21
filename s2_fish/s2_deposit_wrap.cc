#include "deposit.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;

void 
rho_to_rwvars_wrapper(Real_scalar_field &rho, object &zdensity,
            object &xmom, object &ymom)
{
    Array_1d<double> zdensity_array = 
        Array_nd_from_PyObject<double>(zdensity.ptr());
    Array_1d<double> xmom_array = 
        Array_nd_from_PyObject<double>(xmom.ptr());
    Array_1d<double> ymom_array = 
        Array_nd_from_PyObject<double>(ymom.ptr());
    rho_to_rwvars(rho,zdensity_array,xmom_array,ymom_array);
}

void
calculate_rwvars_wrapper(Macro_bunch_store& mbs,
		object &zdensity, object &xmom, object &ymom,
        double z_left, double z_length, object &slice_partition )
{
    Array_1d<double> zdensity_array = 
        Array_nd_from_PyObject<double>(zdensity.ptr());
    Array_1d<double> xmom_array = 
        Array_nd_from_PyObject<double>(xmom.ptr());
    Array_1d<double> ymom_array = 
        Array_nd_from_PyObject<double>(ymom.ptr());
    Array_1d<int> slice_array   =
        Array_nd_from_PyObject<int>(slice_partition.ptr());

    calculate_rwvars(mbs,zdensity_array,xmom_array,ymom_array,
    		z_left,z_length, slice_array);
}

BOOST_PYTHON_MODULE(s2_deposit)
{
    def("deposit_charge_cic", deposit_charge_cic);
    def("deposit_charge_ngp", deposit_charge_ngp);
    def("rho_to_rwvars", rho_to_rwvars_wrapper);
    def("calculate_rwvars", calculate_rwvars_wrapper);
}

