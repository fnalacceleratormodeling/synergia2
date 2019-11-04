

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/simulation/propagator.h"

namespace py = pybind11;
using namespace py::literals;


namespace pybind11
{
    namespace detail
    {
        template <> struct type_caster<std::list<Lattice_element>>
        {
            using value_conv = make_caster<Lattice_element>;

            bool load(handle src, bool convert) 
            { return false; }

            template <typename T>
            static handle cast(T && src, 
                    return_value_policy policy, 
                    handle parent) 
            {
                if (!std::is_lvalue_reference<T>::value)
                    policy = return_value_policy_override<Lattice_element>::policy(policy);
                tuple t(src.size());
                size_t index = 0;

                for (auto &&value : src) 
                {
                    auto value_ = reinterpret_steal<object>(
                            value_conv::cast(forward_like<T>(value), policy, parent));

                    if (!value_)
                        return handle();

                    // steals a reference
                    PyTuple_SET_ITEM(t.ptr(), (ssize_t) index++, value_.release().ptr());
                }

                return t.release();
            }

            PYBIND11_TYPE_CASTER(std::list<Lattice_element>, _("Tuple[") + value_conv::name + _("]"));
        };
    }
}



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

    py::class_<Lattice>(m, "Lattice")
        .def( py::init<>(),
                "Construct an unnamed empty lattice" )
        .def( py::init<std::string const&>(),
                "Construct an empty latttice",
                "name"_a )
        .def( "append",
                &Lattice::append,
                "Append a lattice element to the lattice",
                "element"_a )
        .def( "get_elements",
                py::overload_cast<>(&Lattice::get_elements),
                py::return_value_policy::reference_internal,
                "Get the list of all lattice elements" )
        .def( "get_elements_const",
                py::overload_cast<>(&Lattice::get_elements, py::const_),
                py::return_value_policy::reference_internal,
                "Get the list of all lattice elements" )
        .def( "__repr__", &Lattice::as_string )
        ;
}


