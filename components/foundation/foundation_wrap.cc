#include "four_momentum.h"
#include "reference_particle.h"
#include "distribution.h"
#include <boost/python.hpp>
#include "utils/numpy_multi_ref_converter.h"
#include "utils/comm_converter.h"

using namespace boost::python;

void
(Random_distribution::*random_fill_uniform_ref)(MArray1d_ref, double,
        double) = &Random_distribution::fill_uniform;
void
(Random_distribution::*random_fill_unit_gaussian_ref)(MArray1d_ref) =
    &Random_distribution::fill_unit_gaussian;
void
(Random_distribution::*random_fill_unit_disk_ref)(MArray1d_ref,
        MArray1d_ref) = &Random_distribution::fill_unit_disk;

BOOST_PYTHON_MODULE(pyfoundation)
{
    import_array();
    if (import_mpi4py() < 0) {
        return;
    }
    comm_converter::register_to_and_from_python();
    numpy_multi_array_ref_converter<double,1 >::register_to_and_from_python();
    numpy_const_multi_array_ref_converter<double,1 >::register_to_and_from_python();

    class_<Four_momentum>("Four_momentum", init<double>())
        .def(init<double, double>())
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
        .def("equal",&Four_momentum::equal)
        ;

    class_<Reference_particle>("Reference_particle", init<double, double>())
        .def(init<Four_momentum const &>())
        .def(init<Four_momentum const &,Const_MArray1d_ref>())
        .def("set_four_momentum",&Reference_particle::set_four_momentum)
        .def("set_state",&Reference_particle::set_state)
        .def("set_total_energy",&Reference_particle::set_total_energy)
        .def("get_four_momentum",&Reference_particle::get_four_momentum,
                return_internal_reference<>())
        .def("get_state",&Reference_particle::get_state)
        .def("get_beta",&Reference_particle::get_beta)
        .def("get_gamma",&Reference_particle::get_gamma)
        .def("get_momentum",&Reference_particle::get_momentum)
        .def("get_total_energy",&Reference_particle::get_total_energy)
        .def("equal",&Reference_particle::equal)
        ;

    class_<Distribution, boost::noncopyable > ("Distribution", no_init);

//    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_default_seed_overloads, get_default_seed, 0, 1);
    scope
    Random_distribution_scope =
        class_<Random_distribution, bases<Distribution > > ("Random_distribution",
                init<unsigned long int, Commxx const &>())
            .def(init<unsigned long int, Commxx const &, Random_distribution::Generator>())
            .def("get_default_seed",&Random_distribution::get_default_seed)
            .staticmethod("get_default_seed")
            .def("set_seed",&Random_distribution::set_seed)
            .def("get_original_seed",&Random_distribution::get_original_seed)
//            .def("get_original_seed",(void(*)(const char*))0,get_default_seed_overloads())
            .def("get",&Random_distribution::get)
            .def("fill_uniform",random_fill_uniform_ref)
            .def("fill_unit_gaussian",random_fill_unit_gaussian_ref)
            .def("fill_unit_disk",random_fill_unit_disk_ref)
            ;

        enum_<Random_distribution::Generator >("Generator")
            .value("ranlxd2", Random_distribution::ranlxd2)
            .value("mt19937", Random_distribution::mt19937)
            .export_values();
}
