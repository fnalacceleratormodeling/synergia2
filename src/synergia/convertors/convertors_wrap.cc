#include <boost/python.hpp>
#include "boost/multi_array.hpp"
#include "utils/multi_array_typedefs.h"
#include "utils/numpy_multi_ref_converter.h"
#include "utils/comm_converter.h"
#include "utils/multi_array_serialization.h"
#include "utils/xml_serialization.h"
#include "utils/container_conversions.h"

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

    def("xml_save_array1d", xml_save<MArray1d_ref > );
    def("xml_save_array2d", xml_save<MArray2d_ref > );

    to_python_converter<std::vector<double >,
             container_conversions::to_tuple<std::vector<double > > >();

    to_python_converter<std::list<std::string >,
             container_conversions::to_tuple<std::list<std::string > > >();

}
