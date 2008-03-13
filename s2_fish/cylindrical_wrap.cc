#include "cylindrical.h"
#include "field_domain.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "array_nd/array_nd_python.h"
#include "container_conversions.h"

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
    //---------------------------------------------------------------------
    // std::vector<> conversions
    //---------------------------------------------------------------------
    //~ scitbx::boost_python::container_conversions::from_python_sequence <
    //~ std::vector<int>,
    //~ scitbx::boost_python::container_conversions::variable_capacity_policy > ();

    //~ boost::python::to_python_converter <
    //~ std::vector<int>,
    //~ scitbx::boost_python::container_conversions::to_tuple <
    //~ std::vector<int> > > ();

    //~ scitbx::boost_python::container_conversions::from_python_sequence <
    //~ std::vector<double>,
    //~ scitbx::boost_python::container_conversions::variable_capacity_policy > ();

    //~ boost::python::to_python_converter <
    //~ std::vector<double>,
    //~ scitbx::boost_python::container_conversions::to_tuple <
    //~ std::vector<double> > > ();

    scitbx::boost_python::container_conversions::from_python_sequence <
    std::vector<bool>,
    scitbx::boost_python::container_conversions::variable_capacity_policy > ();

    boost::python::to_python_converter <
    std::vector<bool>,
    scitbx::boost_python::container_conversions::to_tuple <
    std::vector<bool> > > ();

    class_<Field_domain>("Field_domain", init<>())
    .def(init<const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<int> &,
        const std::vector<bool> &>())
    .def("set_params",&Field_domain::set_params)
    .def("get_grid_shape",&Field_domain::get_grid_shape)
    .def("get_cell_size",&Field_domain::get_cell_size)
    .def("get_periodic",&Field_domain::get_periodic)
        //~ .def("get_leftmost_indices_offsets",&Field_domain::get_leftmost_indices_offsets)
    ;
        
    def("get_cylindrical_coords",get_cylindrical_coords_wrapper);
}

