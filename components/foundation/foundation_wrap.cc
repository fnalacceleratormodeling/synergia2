#include "four_momentum.h"
#include "reference_particle.h"
#include <boost/python.hpp>
#include "utils/numpy_multi_ref_converter.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(pyfoundation)
{
    import_array();
    numpy_multi_array_ref_converter<double,1 >::register_to_and_from_python();
    numpy_const_multi_array_ref_converter<double,1 >::register_to_and_from_python();

    class_<Four_momentum>("Four_momentum", init<double>())
        .def("set_total_energy",&Four_momentum::set_total_energy)
        .def("set_kinetic_energy",&Four_momentum::set_kinetic_energy)
        .def("set_momentum",&Four_momentum::set_momentum)
        .def("set_gamma",&Four_momentum::set_gamma)
        .def("set_beta",&Four_momentum::set_beta)
        .def("get_mass",&Four_momentum::get_mass)
        .def("get_total_energy",&Four_momentum::get_total_energy)
        .def("get_kinetic_energy",&Four_momentum::get_kinetic_energy)
        .def("get_momentum",&Four_momentum::get_momentum)
        .def("get_gamma",&Four_momentum::get_gamma)
        .def("get_beta",&Four_momentum::get_beta)
        ;

    class_<Reference_particle>("Reference_particle",
            init<Four_momentum const &>())
        .def(init<Four_momentum const &,Const_MArray1d_ref>())
        .def("set_four_momentum",&Reference_particle::set_four_momentum)
        .def("set_state",&Reference_particle::set_state)
        .def("get_four_momentum",&Reference_particle::get_four_momentum,
                return_internal_reference<>())
        .def("get_state",&Reference_particle::get_state)
        ;
}
