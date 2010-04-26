#include "lattice_element.h"
#include "element_adaptor.h"
#include "lattice.h"
#include <boost/python.hpp>
#include "utils/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(pylattice)
{
    class_<Lattice_element, Lattice_element_sptr >("Lattice_element",
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

    to_python_converter<std::list<Lattice_element_sptr >,
             container_conversions::to_tuple<std::list<Lattice_element_sptr > > >();

    class_<Element_adaptor, Element_adaptor_sptr >("Element_adaptor", init<>())
            .def("set_double_default", &Element_adaptor::set_double_default)
            .def("set_string_default", &Element_adaptor::set_string_default)
            .def("set_default_attributes", &Element_adaptor::set_default_attributes)
            .def("get_chef_elements", &Element_adaptor::get_chef_elements)
            ;

    class_<Element_adaptor_map >("Element_adaptor_map", init<>())
            .def("set_adaptor", &Element_adaptor_map::set_adaptor)
            .def("has_adaptor", &Element_adaptor_map::has_adaptor)
            .def("get_adaptor", &Element_adaptor_map::get_adaptor,
                    return_value_policy<copy_non_const_reference >())
            .def("get_adaptor_names", &Element_adaptor_map::get_adaptor_names)
            ;



    class_<Lattice >("Lattice", init<std::string const& >())
            .def(init<std::string const&, Element_adaptor_map_sptr >())
            .def("get_name", &Lattice::get_name,
                    return_value_policy<copy_const_reference>())
            .def("get_element_adaptor_map", &Lattice::get_element_adaptor_map,
                    return_value_policy<return_by_value >()) // jfa: not sure return_by_value is correct
            .def("set_reference_particle", &Lattice::set_reference_particle)
            .def("has_reference_particle", &Lattice::has_reference_particle)
            .def("get_reference_particle", &Lattice::get_reference_particle,
                    return_value_policy<return_by_value >()) // jfa: not sure return_by_value is correct
            .def("append", &Lattice::append)
            .def("get_elements", &Lattice::get_elements,
                    return_value_policy<copy_non_const_reference >())
            .def("get_length", &Lattice::get_length)
            .def("get_total_angle", &Lattice::get_total_angle)
            .def("print_", &Lattice::print)
            ;
}

