

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/utils/container_conversions.h"

#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"

namespace py = pybind11;
using namespace py::literals;

// convert the std::list<Lattice_element> to a python tuple object
template <> 
struct pybind11::detail::type_caster<std::list<Lattice_element>> 
    : pybind11::detail::py_tuple_caster<std::list<Lattice_element>, Lattice_element> { };


// lattice python module
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
        .value("hkicker",    element_type::hkicker)
        .value("vkicker",    element_type::vkicker)
        .value("sextupole",  element_type::sextupole)
        .value("octupole",   element_type::octupole)
        .value("monitor",    element_type::monitor)
        .value("marker",     element_type::marker)
        .value("instrument", element_type::instrument)
        .value("rcollimator",element_type::rcollimator)
        ;

    py::enum_<marker_type>(m, "marker_type", py::arithmetic())
        .value("h_tunes_corrector", marker_type::h_tunes_corrector)
        .value("v_tunes_corrector", marker_type::v_tunes_corrector)
        .value("h_chrom_corrector", marker_type::h_chrom_corrector)
        .value("v_chrom_corrector", marker_type::v_chrom_corrector)
        ;

    // Lattice_element
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

        .def( "get_length",
                &Lattice_element::get_length )

        // double attribute
        .def( "has_double_attribute", 
                &Lattice_element::has_double_attribute,
                "Check for existence of the named double attribute",
                "name"_a )

        .def( "get_double_attribute", 
                (double (Lattice_element::*)(std::string const&) const)&Lattice_element::get_double_attribute,
                //py::overload_cast<std::string const&>(&Lattice_element::get_double_attribute, py::const_),
                "Get the value of the named double attribute",
                "name"_a )

        .def( "get_double_attribute", 
                (double (Lattice_element::*)(std::string const&, double) const)&Lattice_element::get_double_attribute,
                //py::overload_cast<std::string const&, double>(&Lattice_element::get_double_attribute, py::const_),
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
                (std::string const& (Lattice_element::*)(std::string const&) const)&Lattice_element::get_string_attribute,
                //py::overload_cast<std::string const&>(&Lattice_element::get_string_attribute, py::const_),
                "Get the value of the named double attribute",
                "name"_a )

        .def( "get_string_attribute", 
                (std::string const& (Lattice_element::*)(std::string const&, std::string const&) const)&Lattice_element::get_string_attribute,
                //py::overload_cast<std::string const&, std::string const&>(&Lattice_element::get_string_attribute, py::const_),
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
                (std::vector<double> const& (Lattice_element::*)(std::string const&) const)&Lattice_element::get_vector_attribute,
                //py::overload_cast<std::string const&>(&Lattice_element::get_vector_attribute, py::const_),
                "Get the value of the named vector attribute",
                "name"_a )

        .def( "get_vector_attribute", 
                (std::vector<double> const& (Lattice_element::*)(std::string const&, std::vector<double> const&) const)&Lattice_element::get_vector_attribute,
                //py::overload_cast<std::string const&, std::vector<double> const&>(&Lattice_element::get_vector_attribute, py::const_),
                "Get the value of the named vector attribute, or return the default",
                "name"_a, "default"_a )

        .def( "set_vector_attribute", 
                &Lattice_element::set_vector_attribute, 
                "Set the value of the named vector attribute",
                "name"_a, "value"_a, "increment_revision"_a = true )

        .def( "get_string_attributes",
                &Lattice_element::get_string_attributes )

        .def( "get_double_attributes",
                &Lattice_element::get_double_attributes )

        .def( "set_marker",
                &Lattice_element::set_marker,
                "marker"_a )

        .def( "reset_marker",
                &Lattice_element::reset_marker,
                "marker"_a )

        .def( "has_marker",
                &Lattice_element::has_marker,
                "marker"_a )

        .def( "add_ancestor",
                &Lattice_element::add_ancestor,
                "ancestor"_a )

        .def_static( "get_all_type_names",
                &Lattice_element::get_all_type_names )

        // print
        .def( "print_", 
                &Lattice_element::print,
                "Print the lattice element" )

        .def( "__repr__", &Lattice_element::as_string )
        ;

    // Lattice_element_slice
    py::class_<Lattice_element_slice>(m, "Lattice_element_slice")
        .def( "is_whole", 
                &Lattice_element_slice::is_whole, 
                "Is a whole element" )

        .def( "has_left_edge", 
                &Lattice_element_slice::has_left_edge, 
                "Does this slice include the left edge of the element" )

        .def( "has_right_edge", 
                &Lattice_element_slice::has_right_edge, 
                "Does this slice include the right edge of the element" )

        .def( "get_left", 
                &Lattice_element_slice::get_left, 
                "Get the start position of the slice" )

        .def( "get_right", 
                &Lattice_element_slice::get_right, 
                "Get the end position of the slice" )
        ;

    // Lattice
    py::class_<Lattice>(m, "Lattice")
        .def( py::init<>(), "Construct an unnamed empty lattice" )
        .def( py::init<std::string const&>(), "Construct an empty latttice", "name"_a )
        .def( py::init<Lsexpr const&>(), "Construct from the Lsexpr representation", "lsexpr"_a )

        .def( "get_name",
                &Lattice::get_name,
                "Get the lattice name" )

        .def( "append",
                &Lattice::append,
                "Append a lattice element to the lattice",
                "element"_a )

        .def( "get_reference_particle",
                (Reference_particle& (Lattice::*)())&Lattice::get_reference_particle,
                py::return_value_policy::reference_internal,
                "Get the lattice reference particle" )

        .def( "set_reference_particle",
                &Lattice::set_reference_particle )

        .def( "get_elements",
                (std::list<Lattice_element>& (Lattice::*)())&Lattice::get_elements,
                //py::overload_cast<>(&Lattice::get_elements),
                py::return_value_policy::reference_internal,
                "Get the list of all lattice elements" )

        .def( "get_length",
                &Lattice::get_length,
                "Get the combined length of all elements in the lattice" )

        .def( "get_total_angle",
                &Lattice::get_total_angle,
                "Get the total angle in radians subtended by all elements in the lattice" )

        .def( "get_elements_const",
                (std::list<Lattice_element> const& (Lattice::*)() const)&Lattice::get_elements,
                //py::overload_cast<>(&Lattice::get_elements, py::const_),
                py::return_value_policy::reference_internal,
                "Get the list of all lattice elements" )

        .def( "set_all_double_attribute",
                &Lattice::set_all_double_attribute,
                "Set the value of the named double attribute on all elements",
                "name"_a, "value"_a, "increment_revision"_a = true )

        .def( "set_all_string_attribute",
                &Lattice::set_all_string_attribute,
                "Set the value of the named string attribute on all elements",
                "name"_a, "value"_a, "increment_revision"_a = true )

        .def( "as_json",
                &Lattice::as_json,
                "Returns the lattice object in a json string" )

        .def_static( "load_from_json",
                &Lattice::load_from_json,
                "Create a lattice object from the json string",
                "json_str"_a )

        .def( "__repr__", &Lattice::as_string )
        ;


    using madx_reader_get_lattice_1 = 
        Lattice (MadX_reader::*)(std::string const&);

    using madx_reader_get_lattice_2 = 
        Lattice (MadX_reader::*)(std::string const&, std::string const&);

    // MadX_reader
    py::class_<MadX_reader>(m, "MadX_reader")
        .def( py::init<>(), "Construct a MadX_reader" )
        .def( "get_lattice",
                (madx_reader_get_lattice_1)&MadX_reader::get_lattice,
                "Get the named lattice from an already parsed lattice",
                "line_name"_a )
        .def( "get_lattice",
                (madx_reader_get_lattice_2)&MadX_reader::get_lattice,
                "Parse and get the named lattice",
                "line_name"_a,
                "filename"_a )
        .def( "parse",
                &MadX_reader::parse,
                "Parse a lattice string",
                "lattice"_a )
        .def( "parse_file",
                &MadX_reader::parse_file,
                "Parse a lattice file",
                "filename"_a )
        ;
}


