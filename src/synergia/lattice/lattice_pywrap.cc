

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/simulation/propagator.h"

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(lattice, m)
{
    py::enum_<element_type>(m, "element_type", py::arithmetic())
        .value("generic",    element_type::generic)
        .value("drift",      element_type::drift)
        .value("rbend",      element_type::rbend)
        .value("sbend",      element_type::sbend)
        .value("quadrupole", element_type::quadrupole)
        .value("multipole",  element_type::multipole)
        .value("rfcavity",   element_type::rfcavity)
        ;

    py::class_<Lattice_element>(m, "Lattice_element")
        .def( py::init<>(),
                "Construct a generic lattice element" )

        .def( py::init<std::string const&, std::string const&>(), 
                "Construct a lattice element with type and name",
                "type"_a, "name"_a )

        .def( "get_type", 
                &Lattice_element::get_type, 
                "Returns lattice element type" )

        .def( "get_type_name", 
                &Lattice_element::get_type_name, 
                "Returns lattice element type string" )

        .def( "get_name", 
                &Lattice_element::get_name, 
                "Returns lattice element name" )

        // double attribute
        .def( "has_double_attribute", 
                &Lattice_element::has_double_attribute,
                "Check for existence of the named double attribute",
                "name"_a )

        .def( "get_double_attribute", 
                py::overload_cast<std::string const&>(&Lattice_element::get_double_attribute, py::const_),
                "Get the value of the named double attribute",
                "name"_a )

        .def( "get_double_attribute", 
                py::overload_cast<std::string const&, double>(&Lattice_element::get_double_attribute, py::const_),
                "Get the value of the named double attribute, or return the default",
                "name"_a, "default"_a )

        .def( "set_double_attribute", 
                &Lattice_element::set_double_attribute, 
                "Set the value of the named double attribute", 
                "name"_a, "value"_a, "increment_revision"_a = true )

        // string attribute
        .def( "has_string_attribute", 
                &Lattice_element::has_string_attribute,
                "Check for existence of the named string attribute",
                "name"_a )

        .def( "get_string_attribute", 
                py::overload_cast<std::string const&>(&Lattice_element::get_string_attribute, py::const_),
                "Get the value of the named double attribute",
                "name"_a )

        .def( "get_string_attribute", 
                py::overload_cast<std::string const&, std::string const&>(&Lattice_element::get_string_attribute, py::const_),
                "Get the value of the named string attribute, or return the default",
                "name"_a, "default"_a )

        .def( "set_string_attribute", 
                &Lattice_element::set_string_attribute, 
                "Set the value of the named string attribute", 
                "name"_a, "value"_a, "increment_revision"_a = true )

        // vector attribute
        .def( "has_vector_attribute", 
                &Lattice_element::has_vector_attribute,
                "Check for existence of the named vector attribute",
                "name"_a )

        .def( "get_vector_attribute", 
                py::overload_cast<std::string const&>(&Lattice_element::get_vector_attribute, py::const_),
                "Get the value of the named vector attribute",
                "name"_a )

        .def( "get_vector_attribute", 
                py::overload_cast<std::string const&, std::vector<double> const&>(&Lattice_element::get_vector_attribute, py::const_),
                "Get the value of the named vector attribute, or return the default",
                "name"_a, "default"_a )

        .def( "set_vector_attribute", 
                &Lattice_element::set_vector_attribute, 
                "Set the value of the named vector attribute",
                "name"_a, "value"_a, "increment_revision"_a = true )

        // print
        .def( "print_", 
                &Lattice_element::print,
                "Print the lattice element" )

        .def( "__repr__", &Lattice_element::as_string )
        ;
}


