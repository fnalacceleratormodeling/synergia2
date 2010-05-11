#include "four_momentum.h"
#include "reference_particle.h"
#include "distribution.h"
#include <boost/python.hpp>
#include "utils/numpy_multi_ref_converter.h"
#include "utils/comm_converter.h"
#include "math_constants.h"
#include "physical_constants.h"

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

class Dummy
{

};

BOOST_PYTHON_MODULE(pyfoundation)
{
    import_array();
    if (import_mpi4py() < 0) {
        return;
    }

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

    class_<Reference_particle>("Reference_particle",
            init<int, double, double>())
        .def(init<int, Four_momentum const &>())
        .def(init<int, Four_momentum const &,Const_MArray1d_ref>())
        .def("set_four_momentum",&Reference_particle::set_four_momentum)
        .def("set_state",&Reference_particle::set_state)
        .def("set_total_energy",&Reference_particle::set_total_energy)
        .def("increment_trajectory", &Reference_particle::increment_trajectory)
        .def("start_repetition", &Reference_particle::start_repetition)
        .def("set_trajectory", &Reference_particle::set_trajectory)
        .def("get_charge", &Reference_particle::get_charge)
        .def("get_four_momentum",&Reference_particle::get_four_momentum,
                return_internal_reference<>())
        .def("get_state",&Reference_particle::get_state)
        .def("get_beta",&Reference_particle::get_beta)
        .def("get_gamma",&Reference_particle::get_gamma)
        .def("get_momentum",&Reference_particle::get_momentum)
        .def("get_total_energy",&Reference_particle::get_total_energy)
        .def("get_trajectory_length", &Reference_particle::get_trajectory_length)
        .def("get_s", &Reference_particle::get_s)
        .def("get_repetition", &Reference_particle::get_repetition)
        .def("get_repetition_length", &Reference_particle::get_repetition_length)
        .def("equal",&Reference_particle::equal)
        ;

    class_<Distribution, boost::noncopyable > ("Distribution", no_init);

    class_<Dummy >("constants",no_init)
    // math_constants.h
        .def_readonly("pi", constants::pi)
    // physical_constants.h
        .def_readonly("mp", constants::mp)
        .def_readonly("me", constants::me)
        .def_readonly("mmu", constants::mmu)
        .def_readonly("proton_charge", constants::proton_charge)
        .def_readonly("antiproton_charge", constants::antiproton_charge)
        .def_readonly("electron_charge", constants::electron_charge)
        .def_readonly("positron_charge", constants::positron_charge)
        .def_readonly("muon_charge", constants::muon_charge)
        .def_readonly("antimuon_charge", constants::antimuon_charge)
        ;

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
