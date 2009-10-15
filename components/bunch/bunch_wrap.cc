#include "bunch.h"
#include <boost/python.hpp>
#include "utils/numpy_multi_ref_converter.h"
#include "utils/comm_converter.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(pybunch)
{
    import_array();
    if (import_mpi4py() < 0) {
        return;
    }
    comm_converter::register_to_and_from_python();
    numpy_multi_array_ref_converter<double, 2 >::register_to_and_from_python();

    class_<Fixed_t_z_converter, boost::noncopyable > ("Fixed_t_z_converter",
            no_init);
    class_<Fixed_t_z_zeroth, bases<Fixed_t_z_converter > > ("Fixed_t_z_zeroth",
            init< > ());
    class_<Fixed_t_z_ballistic, bases<Fixed_t_z_converter > > (
            "Fixed_t_z_ballistic", init< > ());

    scope
            Bunch_scope =
                    class_<Bunch > ("Bunch", init<Reference_particle const&,
                            int, int, double, Commxx const& > ())
                            .def("set_particle_charge", &Bunch::set_particle_charge)
                            .def("set_real_num", &Bunch::set_real_num)
                            .def("set_local_num", &Bunch::set_local_num)
                            .def("update_total_num", &Bunch::update_total_num)
                            .def("set_converter", &Bunch::set_converter)
                            .def("convert_to_state", &Bunch::convert_to_state)
                            .def("get_reference_particle",
                                    &Bunch::get_reference_particle,
                                    return_internal_reference< > ())
                            .def("get_local_particles", &Bunch::get_local_particles)
                            .def("get_particle_charge", &Bunch::get_particle_charge)
                            .def("get_mass", &Bunch::get_mass)
                            .def("get_real_num", &Bunch::get_real_num)
                            .def("get_local_num", &Bunch::get_local_num)
                            .def("get_total_num", &Bunch::get_total_num)
                            .def("get_state", &Bunch::get_state)
                            .def("inject", &Bunch::inject)
                            ;

    enum_<Bunch::State > ("State") .value("fixed_z", Bunch::fixed_z) .value(
            "fixed_t", Bunch::fixed_t) .export_values();

    // Ideally, this block would use, e.g., Bunch::x instead of 0. Somehow,
    // however, that doesn't work.
    scope().attr("x") = 0;
    scope().attr("xp") = 1;
    scope().attr("y") = 2;
    scope().attr("yp") = 3;
    scope().attr("z") = 4;
    scope().attr("zp") = 5;
    scope().attr("t") = 4;
    scope().attr("tp") = 5;
    scope().attr("id") = 6;
}
