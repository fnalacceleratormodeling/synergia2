#include "electric_field.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(s2_electric_field)
{
    def("calculate_E_n", calculate_E_n);
    def("apply_E_n_kick", apply_E_n_kick);
    def("full_kick", full_kick);
}

