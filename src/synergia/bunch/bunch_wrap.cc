#include "bunch.h"
#include "diagnostics.h"
#include "diagnostics_writer.h"
#include "populate.h"
#include <boost/python.hpp>
#include "synergia/utils/numpy_multi_ref_converter.h"
#include "synergia/utils/comm_converter.h"

using namespace boost::python;

PyObject *
bunch_get_comm_workaround(Bunch const& bunch)
{
    PyObject *retval;
    retval = PyMPIComm_New(bunch.get_comm().get());
    return retval;
}

BOOST_PYTHON_MODULE(bunch)
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

    class_<Diagnostics, Diagnostics_sptr >("Diagnostics",init< >())
        .def(init<Bunch const & >())
        .def("update", &Diagnostics::update)
        .def("get_s", &Diagnostics::get_s)
        .def("get_repetition", &Diagnostics::get_repetition)
        .def("get_trajectory_length", &Diagnostics::get_trajectory_length)
        .def("get_mean", &Diagnostics::get_mean)
        .def("get_std", &Diagnostics::get_std)
        ;

    class_<Diagnostics_full2, Diagnostics_full2_sptr, bases<Diagnostics > >
            ("Diagnostics_full2",init< >())
        .def(init<Bunch const & >())
        .def("get_mom2", &Diagnostics_full2::get_mom2)
        .def("get_corr", &Diagnostics_full2::get_corr)
        .def("get_emitx", &Diagnostics_full2::get_emitx)
        .def("get_emity", &Diagnostics_full2::get_emity)
        .def("get_emitz", &Diagnostics_full2::get_emitz)
        .def("get_emitxy", &Diagnostics_full2::get_emitxy)
        .def("get_emitxyz",&Diagnostics_full2::get_emitxyz)
        ;

    class_<Diagnostics_track, Diagnostics_track_sptr, bases<Diagnostics > >
        ("Diagnostics_track",init<int >())
        .def(init<Bunch const &, int >())
        .def("update", &Diagnostics_track::update)
        ;

    class_<Diagnostics_particles, Diagnostics_particles_sptr, bases<Diagnostics > >
            ("Diagnostics_particles",init< >())
        .def(init<int >())
        .def(init<Bunch const & >())
        .def(init<Bunch const &, int>())
        ;

    class_<Diagnostics_writer, Diagnostics_writer_sptr >("Diagnostics_writer",
            init<std::string const& , Diagnostics_sptr >())
        .def(init< >())
        .def("is_dummy", &Diagnostics_writer::is_dummy)
        .def("get_diagnostics", &Diagnostics_writer::get_diagnostics_sptr)
        .def("write", &Diagnostics_writer::write)
        .def("update_and_write", &Diagnostics_writer::update_and_write)
        ;

    class_<Multi_diagnostics_writer >("Multi_diagnostics_writer",
            init<>())
        .def("append", &Multi_diagnostics_writer::append)
        ;

    def("no_diagnostics", no_diagnostics);
    def("populate_6d", populate_6d);

    typedef Reference_particle & (Bunch::*get_reference_particle_non_const_type)();
    typedef MArray2d_ref (Bunch::*get_local_particles_non_const_type)();
    scope
        Bunch_scope =
            class_<Bunch > ("Bunch", init<Reference_particle const&,
                    int, double, Commxx const& > ())
                .def(init<Reference_particle const&, int, double,
                        Commxx const&, int >())
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
                // jfa: the following implementation does not work for reasons I do
                //      not understand.
//                .def("get_comm",
//                            &Bunch::get_comm, return_internal_reference< > ())
                .def("get_comm",
                            &bunch_get_comm_workaround)
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
