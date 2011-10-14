#include "lattice_element.h"
#include "element_adaptor.h"
#include "lattice.h"
#include "lattice_diagnostics.h"
#include "chef_lattice.h"
#include "chef_utils.h"
#include <boost/python.hpp>
#include <boost/python/dict.hpp>
#include "synergia/utils/container_conversions.h"
#include "synergia/utils/serialization.h"

using namespace boost::python;

// jfa: ideally, this would be done through a container conversion.
// This implementation is simpler, if less general.
dict
get_double_attributes_workaround(Lattice_element const& lattice_element)
{
    dict retval;
    for (std::map<std::string, double >::const_iterator it =
            lattice_element.get_double_attributes().begin(); it
            != lattice_element.get_double_attributes().end(); ++it) {
        retval[it->first] = it->second;
    }
    return retval;
}

// jfa: ideally, this would be done through a container conversion.
// This implementation is simpler, if less general.
dict
get_string_attributes_workaround(Lattice_element const& lattice_element)
{
    dict retval;
    for (std::map<std::string, string >::const_iterator it =
            lattice_element.get_string_attributes().begin(); it
            != lattice_element.get_string_attributes().end(); ++it) {
        retval[it->first] = it->second;
    }
    return retval;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Lattice_element_set_double_attribute_overloads,
        set_double_attribute, 2, 3);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Lattice_element_set_string_attribute_overloads,
        set_string_attribute, 2, 3);

BOOST_PYTHON_MODULE(lattice)
{
//    import("pyconvertors");
    class_<Lattice_element, Lattice_element_sptr >("Lattice_element",
            init<std::string, std::string >())
        .def("get_type", &Lattice_element::get_type,
                return_value_policy<copy_const_reference>())
        .def("get_name", &Lattice_element::get_name,
                return_value_policy<copy_const_reference>())
        .def("add_ancestor", &Lattice_element::add_ancestor)
        .def("get_ancestors", &Lattice_element::get_ancestors,
                return_value_policy<copy_const_reference>())
        .def("set_double_attribute", &Lattice_element::set_double_attribute,
                Lattice_element_set_double_attribute_overloads())
        .def("has_double_attribute", &Lattice_element::has_double_attribute)
        .def("get_double_attribute", &Lattice_element::get_double_attribute)
        .def("get_double_attributes", get_double_attributes_workaround)
        .def("set_string_attribute", &Lattice_element::set_string_attribute,
                Lattice_element_set_string_attribute_overloads())
        .def("has_string_attribute", &Lattice_element::has_string_attribute)
        .def("get_string_attribute", &Lattice_element::get_string_attribute,
                return_value_policy<copy_const_reference>())
        .def("get_string_attributes", get_string_attributes_workaround)
        .def("set_length_attribute_name", &Lattice_element::set_length_attribute_name)
        .def("set_bend_angle_attribute_name", &Lattice_element::set_bend_angle_attribute_name)
        .def("get_length", &Lattice_element::get_length)
        .def("get_bend_angle", &Lattice_element::get_bend_angle)
        .def("get_revision", &Lattice_element::get_revision)
        .def("print_",&Lattice_element::print)
       ;

    to_python_converter<std::list<Lattice_element_sptr >,
             container_conversions::to_tuple<Lattice_elements > >();
    container_conversions::from_python_sequence<Lattice_elements,
             container_conversions::variable_capacity_policy >();

    class_<Element_adaptor, Element_adaptor_sptr >("Element_adaptor", init<>())
            .def("set_double_default", &Element_adaptor::set_double_default)
            .def("set_string_default", &Element_adaptor::set_string_default)
            .def("set_default_attributes", &Element_adaptor::set_default_attributes)
            .def("get_chef_elements", &Element_adaptor::get_chef_elements)
            ;

    class_<Element_adaptor_map >("Element_adaptor_map", init<>())
            .def("set_adaptor", &Element_adaptor_map::set_adaptor)
            .def("has_adaptor", &Element_adaptor_map::has_adaptor)
            .def("get_adaptor", &Element_adaptor_map::get_adaptor)
            .def("get_adaptor_names", &Element_adaptor_map::get_adaptor_names)
            ;

    class_<Lattice, Lattice_sptr >("Lattice", init<std::string const& >())
            .def(init<>())
            .def("get_name", &Lattice::get_name,
                    return_value_policy<copy_const_reference>())
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

    class_<Chef_lattice, Chef_lattice_sptr >("Chef_lattice", init<Lattice_sptr>())
            .def(init<Lattice_sptr, Element_adaptor_map_sptr>())
//            .def("get_chef_elements", &Chef_lattice::get_chef_elements,
//                    return_value_policy<copy_non_const_reference >())
            .def("get_beamline", &Chef_lattice::get_beamline_sptr)
            .def("get_sliced_beamline", &Chef_lattice::get_sliced_beamline_sptr)
            ;

    class_<Lattice_diagnostics, Lattice_diagnostics_sptr, bases<Generalized_diagnostics > >
        ("Lattice_diagnostics",init<Lattice_sptr, std::string const&,
                std::string const& >())
        .def("set_default_value", &Lattice_diagnostics::set_default_value)
        .def("get_default_value", &Lattice_diagnostics::get_default_value)
        .def("set_reduce", &Lattice_diagnostics::set_reduce)
        .def("get_reduce", &Lattice_diagnostics::get_reduce)
//        .def("set_reduce_op", &Lattice_diagnostics::set_reduce_op)
//        .def("get_reduce_op", &Lattice_diagnostics::get_reduce_op)
        .def("update", &Lattice_diagnostics::update)
        .def("write", &Lattice_diagnostics::write)
        .def("update_and_write", &Lattice_diagnostics::update_and_write)
        ;

    def("print_chef_beamline", print_chef_beamline);
    def("reference_particle_to_chef_particle",
            reference_particle_to_chef_particle);
    def("reference_particle_to_chef_jet_particle",
            reference_particle_to_chef_jet_particle);
//    propagate_reference_particle(Reference_particle const& reference_particle,
//            BmlPtr beamline_sptr);
   def("chef_unit_conversion", chef_unit_conversion);
   def("get_chef_index",get_chef_index);

   def("xml_save_lattice", xml_save<Lattice > );
   def("xml_load_lattice", xml_load<Lattice > );

}

