#include "lattice_element.h"
#include "element_adaptor.h"
#include "lattice.h"
#include "chef_lattice.h"
#include "chef_utils.h"
#include <boost/python.hpp>
#include "utils/container_conversions.h"
#include "utils/xml_serialization.h"

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

   to_python_converter<std::vector<double >,
            container_conversions::to_tuple<std::vector<double > > >();

}

