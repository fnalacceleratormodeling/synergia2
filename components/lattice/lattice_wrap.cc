#include "lattice_element.h"
#include <boost/python.hpp>
#include "utils/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(pylattice)
{
    class_<Lattice_element>("Lattice_element",
            init<std::string, std::string >())
        .def("get_type", &Lattice_element::get_type,
                return_value_policy<copy_const_reference>())
        .def("get_name", &Lattice_element::get_name,
                return_value_policy<copy_const_reference>())
        .def("add_ancestor", &Lattice_element::add_ancestor)
        .def("get_ancestors", &Lattice_element::get_ancestors,
                return_value_policy<copy_const_reference>())
        .def("set_double_attribute", &Lattice_element::set_double_attribute)
        .def("has_double_attribute", &Lattice_element::has_double_attribute)
        .def("get_double_attribute", &Lattice_element::get_double_attribute)
//        .def("get_double_attributes", &Lattice_element::get_double_attributes,
//                return_value_policy<copy_const_reference>())
        .def("set_string_attribute", &Lattice_element::set_string_attribute)
        .def("has_string_attribute", &Lattice_element::has_string_attribute)
        .def("get_string_attribute", &Lattice_element::get_string_attribute,
                return_value_policy<copy_const_reference>())
//        .def("get_string_attributes", &Lattice_element::get_string_attributes,
//                return_value_policy<copy_const_reference>())
        .def("set_length_attribute_name", &Lattice_element::set_length_attribute_name)
        .def("set_bend_angle_attribute_name", &Lattice_element::set_bend_angle_attribute_name)
        .def("get_length", &Lattice_element::get_length)
        .def("get_bend_angle", &Lattice_element::get_bend_angle)
       ;

    to_python_converter<std::list<std::string >,
             container_conversions::to_tuple<std::list<std::string > > >();
}

