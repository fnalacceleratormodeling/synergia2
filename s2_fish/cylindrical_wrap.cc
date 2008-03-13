#include "cylindrical.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"

using namespace boost::python;

#include <iostream>

void 
get_cylindrical_coords_wrapper(Macro_bunch_store &mbs, 
    object &coords)
{
    Array_2d<double> coords_array = 
        Array_nd_from_PyObject<double>(coords.ptr());
    get_cylindrical_coords(mbs,coords_array);
}

BOOST_PYTHON_MODULE(s2_solver_cylindrical)
{
    def("get_cylindrical_coords",get_cylindrical_coords_wrapper);
}

