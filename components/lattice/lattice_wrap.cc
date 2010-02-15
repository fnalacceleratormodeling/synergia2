#include "lattice_element.h"
#include <boost/python.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(pylattice)
{
    class_<Lattice_element>("Lattice_element",
            init<std::string, std::string >())
        .def("get_type",&Lattice_element::get_type,
                return_value_policy<copy_const_reference>())
        .def("get_name",&Lattice_element::get_name,
                return_value_policy<copy_const_reference>())
        .def("set_double_attribute",&Lattice_element::set_double_attribute)
        .def("has_double_attribute",&Lattice_element::has_double_attribute)
        .def("get_double_attribute",&Lattice_element::get_double_attribute)
        ;
}

