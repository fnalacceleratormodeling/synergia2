#include "electric_field.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;


BOOST_PYTHON_MODULE(s2_electric_field)
{
    def("calculate_E_n", calculate_E_n);
    def("apply_E_n_kick", apply_E_n_kick);
    def("full_kick", full_kick);
    def("transverse_kick", transverse_kick);
}

