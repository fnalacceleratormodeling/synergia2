#include "bunch.h"
#include "diagnostics.h"
#include "populate.h"
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

    class_<Fixed_t_z_converter, boost::noncopyable > ("Fixed_t_z_converter",
            no_init);
    class_<Fixed_t_z_zeroth, bases<Fixed_t_z_converter > > ("Fixed_t_z_zeroth",
            init< > ());
    class_<Fixed_t_z_ballistic, bases<Fixed_t_z_converter > > (
            "Fixed_t_z_ballistic", init< > ());

    class_<Diagnostics >("Diagnostics",init< >())
        .def(init<Bunch const &, double  >())
        .def("update", &Diagnostics::update)
        .def("get_s", &Diagnostics::get_s)
        .def("get_mean", &Diagnostics::get_mean)
        .def("get_std", &Diagnostics::get_std)
        ;

    class_<Diagnostics_full2, bases<Diagnostics > >("Diagnostics_full2",init< >())
        .def(init<Bunch const &, double >())
        .def("get_mom2",&Diagnostics_full2::get_mom2)
        .def("get_corr",&Diagnostics_full2::get_corr)
        .def("get_emitx",&Diagnostics_full2::get_emitx)
        .def("get_emity",&Diagnostics_full2::get_emity)
        .def("get_emitz",&Diagnostics_full2::get_emitz)
        .def("get_emitxy",&Diagnostics_full2::get_emitxy)
        .def("get_emitxyz",&Diagnostics_full2::get_emitxyz)
        ;

    def("populate_6d", populate_6d);

    typedef Reference_particle & (Bunch::*get_reference_particle_non_const_type)();
    typedef MArray2d_ref (Bunch::*get_local_particles_non_const_type)();
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
                        get_reference_particle_non_const_type(
                                &Bunch::get_reference_particle),
                        return_internal_reference< >())
                .def("get_local_particles",
                        get_local_particles_non_const_type(
                                &Bunch::get_local_particles))
                .def("get_particle_charge", &Bunch::get_particle_charge)
                .def("get_mass", &Bunch::get_mass)
                .def("get_real_num", &Bunch::get_real_num)
                .def("get_local_num", &Bunch::get_local_num)
                .def("get_total_num", &Bunch::get_total_num)
                .def("get_state", &Bunch::get_state)
//                .def("get_comm",
//                            &Bunch::get_comm, return_internal_reference< > ())
                .def("inject", &Bunch::inject)
                ;

    enum_<Bunch::State > ("State")
        .value("fixed_z", Bunch::fixed_z)
        .value("fixed_t", Bunch::fixed_t)
        .export_values();

    scope().attr("x") = Bunch::x;
    scope().attr("xp") = Bunch::xp;
    scope().attr("y") = Bunch::y;
    scope().attr("yp") = Bunch::yp;
    scope().attr("z") = Bunch::z;
    scope().attr("zp") = Bunch::zp;
    scope().attr("t") = Bunch::t;
    scope().attr("tp") = Bunch::tp;
    scope().attr("id") = Bunch::id;
}
