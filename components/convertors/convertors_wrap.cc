#include <boost/python.hpp>
#include "boost/multi_array.hpp"
#include "utils/numpy_multi_ref_converter.h"
#include "utils/comm_converter.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(pyconvertors)
{
    import_array();
    if (import_mpi4py() < 0) {
        return;
    }
    comm_converter::register_to_and_from_python();
    numpy_multi_array_ref_converter<double, 1 >::register_to_and_from_python();
    numpy_const_multi_array_ref_converter<double, 1 >::register_to_and_from_python();
    numpy_multi_array_ref_converter<double, 2 >::register_to_and_from_python();
    numpy_const_multi_array_ref_converter<double, 2 >::register_to_and_from_python();
}
