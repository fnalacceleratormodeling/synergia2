#include <boost/python.hpp>
#include "bunch.h"
#include "bunch_train.h"
#include "diagnostics.h"
#include "diagnostics_basic.h"
#include "diagnostics_full2.h"
#include "diagnostics_particles.h"
#include "diagnostics_track.h"
#include "diagnostics_bulk_track.h"
#include "diagnostics_reference_particle.h"
#include "populate.h"
#include "analysis.h"
#include "synergia/utils/numpy_multi_ref_converter.h"
#include "synergia/utils/comm_converter.h"
#include "synergia/utils/container_conversions.h"
#include "synergia/foundation/multi_diagnostics.h"

using namespace boost::python;

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
  
    class_<Diagnostics, Diagnostics_sptr, boost::noncopyable >
        ("Diagnostics", no_init)
        .def("get_filename", &Diagnostics::get_filename,
                return_value_policy<copy_const_reference >())
        .def("set_bunch", &Diagnostics::set_bunch_sptr)
        .def("have_bunch", &Diagnostics::have_bunch)
        .def("get_bunch", &Diagnostics::get_bunch,
                return_internal_reference< >())
        .def("is_serial", &Diagnostics::is_serial)
        .def("update", &Diagnostics::update)
        .def("write", &Diagnostics::write)
        .def("update_and_write", &Diagnostics::update_and_write)
        ;

    container_conversions::from_python_sequence<std::list<Diagnostics_sptr >,
            container_conversions::variable_capacity_policy >();

    to_python_converter<std::list<Diagnostics_sptr >,
            container_conversions::to_tuple<std::list<Diagnostics_sptr > > >();

    class_<Diagnostics_basic, Diagnostics_basic_sptr, bases<Diagnostics > >
        ("Diagnostics_basic",init<std::string const& >())
        .def(init<std::string const&, std::string const& >())
        .def("set_bunch", &Diagnostics_basic::set_bunch_sptr)
        .def("get_s_n", &Diagnostics_basic::get_s_n)
        .def("get_repetition", &Diagnostics_basic::get_repetition)
        .def("get_s", &Diagnostics_basic::get_s)
        .def("get_mean", &Diagnostics_basic::get_mean)
        .def("get_std", &Diagnostics_basic::get_std)
        .def("get_min", &Diagnostics_basic::get_min)
        .def("get_max", &Diagnostics_basic::get_max)
        ;

    class_<Diagnostics_full2, Diagnostics_full2_sptr, bases<Diagnostics > >
        ("Diagnostics_full2",init<std::string const& >())
        .def(init<std::string const& , std::string const& >())
        .def("set_bunch", &Diagnostics_full2::set_bunch_sptr)
        .def("get_s_n", &Diagnostics_full2::get_s_n)
        .def("get_repetition", &Diagnostics_full2::get_repetition)
        .def("get_s", &Diagnostics_full2::get_s)
        .def("get_mean", &Diagnostics_full2::get_mean)
        .def("get_std", &Diagnostics_full2::get_std)
        .def("get_min", &Diagnostics_basic::get_min)
        .def("get_max", &Diagnostics_basic::get_max)
        .def("get_mom2", &Diagnostics_full2::get_mom2)
        .def("get_corr", &Diagnostics_full2::get_corr)
        .def("get_emitx", &Diagnostics_full2::get_emitx)
        .def("get_emity", &Diagnostics_full2::get_emity)
        .def("get_emitz", &Diagnostics_full2::get_emitz)
        .def("get_emitxy", &Diagnostics_full2::get_emitxy)
        .def("get_emitxyz",&Diagnostics_full2::get_emitxyz)
        ;

    class_<Diagnostics_track, Diagnostics_track_sptr, bases<Diagnostics > >
    		("Diagnostics_track",init<std::string const&, int>())
    	.def(init<std::string const&, int, std::string const& >())
        .def("set_bunch", &Diagnostics_track::set_bunch_sptr)
        ;

    class_<Diagnostics_bulk_track, Diagnostics_bulk_track_sptr, bases<Diagnostics> >
      ("Diagnostics_bulk_track",init<std::string const&, int>())
      .def(init<std::string const&, int, int >())
      .def(init<std::string const&, int, int, std::string const& >())
      .def("set_bunch", &Diagnostics_bulk_track::set_bunch_sptr)
      ;

    class_<Diagnostics_particles, Diagnostics_particles_sptr, bases<Diagnostics > >
            ("Diagnostics_particles",init<std::string const&>())
        .def(init<std::string const&, int >())
        .def(init<std::string const&, int, int >())
        .def(init<std::string const&, int, int, std::string const& >())
        .def(init<std::string const&, int >())
        .def(init<std::string const&, int, int>())
        .def(init<std::string const&, int, int, std::string const& >())
        .def("set_bunch", &Diagnostics_particles::set_bunch_sptr)
        ;

    class_<Diagnostics_reference_particle, Diagnostics_reference_particle_sptr, bases<Diagnostics > >
        ("Diagnostics_reference_particle",init<std::string const& , std::string const& >())
        .def("set_bunch", &Diagnostics_reference_particle::set_bunch_sptr)
        ;

    class_<Multi_diagnostics>("Multi_diagnostics",
            init<>())
        .def("append", &Multi_diagnostics::append)
        ;

    class_<Core_diagnostics >
    ("Core_diagnostics", init<>())
    		.def("calculate_mean", &Core_diagnostics::calculate_mean)
    		.staticmethod("calculate_mean")
    		.def("calculate_z_mean", &Core_diagnostics::calculate_z_mean)
    		.staticmethod("calculate_z_mean")
    		.def("calculate_spatial_mean", &Core_diagnostics::calculate_spatial_mean)
    		.staticmethod("calculate_spatial_mean")
    		.def("calculate_std", &Core_diagnostics::calculate_std)
    		.staticmethod("calculate_std")
    		.def("calculate_spatial_std", &Core_diagnostics::calculate_spatial_std)
    		.staticmethod("calculate_spatial_std")
    		.def("calculate_mom2", &Core_diagnostics::calculate_mom2)
    		.staticmethod("calculate_mom2")
    		.def("calculate_min", &Core_diagnostics::calculate_min)
    		.staticmethod("calculate_min")
    		.def("calculate_max", &Core_diagnostics::calculate_max)
    		.staticmethod("calculate_max")
            .def("print_bunch_parameters", &Core_diagnostics::print_bunch_parameters)
            .staticmethod("print_bunch_parameters")
    		;



    class_<Analysis>("Analysis", init<std::string const&, int>())
      .def("get_betatron_averaged_tune", &Analysis::get_betatron_averaged_tune)
      .def("get_betatron_tune_spread", &Analysis::get_betatron_tune_spread)
      .def("get_betatron_tunes", &Analysis::get_betatron_tunes)
      .def("get_XCoords_forTunes", &Analysis::get_XCoords_forTunes)
      .def("get_YCoords_forTunes", &Analysis::get_YCoords_forTunes)
      .def("get_num_betatron_tune_found", &Analysis::get_num_betatron_tune_found)
      .def("set_minimum_required_turn_Num", &Analysis::set_minimum_required_turn_Num)
      .def("get_minimum_required_turn_Num", &Analysis::get_minimum_required_turn_Num)
      .def("get_num_turns", &Analysis::get_num_turns)
      .def("get_num_part_first_bunch", &Analysis::get_num_part_first_bunch)
      .def("get_transverse_action_for_particle", &Analysis::get_transverse_action_for_particle)
      .def("get_transverse_action_for_bunch", &Analysis::get_transverse_action_for_bunch)
    ;

    def("no_diagnostics", no_diagnostics);
    def("populate_6d", populate_6d);
    def("populate_6d_truncated", populate_6d_truncated);
    def("populate_transverse_gaussian", populate_transverse_gaussian);
    def("populate_uniform_cylinder", populate_uniform_cylinder);
    def("populate_transverse_KV_GaussLong", populate_transverse_KV_GaussLong);
    def("populate_two_particles", populate_two_particles);
    def("get_correlation_matrix",get_correlation_matrix);
    def("populate_transverse_KV_GaussLong", populate_transverse_KV_GaussLong);
    def("populate_transverseKV_logitudinalGaussian", populate_transverseKV_logitudinalGaussian);
    def("populate_longitudinal_boxcar", populate_longitudinal_boxcar);
    def("populate_longitudinal_uniform",populate_longitudinal_uniform);

    container_conversions::from_python_sequence<std::vector<Bunch_sptr >,
            container_conversions::variable_capacity_policy >();

    to_python_converter<std::vector<Bunch_sptr >,
            container_conversions::to_tuple<std::vector<Bunch_sptr > > >();

    class_<Bunch_train, Bunch_train_sptr >("Bunch_train",
            init<Bunches, double >())
            .def(init<Bunches, std::vector<double > const& >())
            .def("get_size", &Bunch_train::get_size)
            .def("get_bunches", &Bunch_train::get_bunches, return_value_policy<copy_non_const_reference >())
            .def("get_spacings", &Bunch_train::get_spacings, return_value_policy<copy_non_const_reference >())
            ;

    typedef Reference_particle & (Bunch::*get_reference_particle_non_const_type)();
    typedef MArray2d_ref (Bunch::*get_local_particles_non_const_type)();
    scope
        Bunch_scope =
            class_<Bunch, Bunch_sptr > ("Bunch", init<Reference_particle const&,
                    int, double, Commxx_sptr > ())                        
                .def("set_particle_charge", &Bunch::set_particle_charge)
                .def("set_real_num", &Bunch::set_real_num)
                .def("set_local_num", &Bunch::set_local_num)
                .def("set_bucket_index", &Bunch::set_bucket_index)
                .def("set_z_period_length", &Bunch::set_z_period_length)
                .def("set_longitudinal_aperture_length", &Bunch::set_longitudinal_aperture_length)
                .def("set_bucket_barrier_length", &Bunch::set_bucket_barrier_length)
                .def("set_open_longitudinal_boundary", &Bunch::set_open_longitudinal_boundary)
                .def("set_longitudinal_boundary", &Bunch::set_longitudinal_boundary)
                .def("update_total_num", &Bunch::update_total_num)
                .def("set_sort_period", &Bunch::set_sort_period)
                .def("sort", &Bunch::sort)
                .def("periodic_sort", &Bunch::periodic_sort)
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
                .def("get_z_period_length", &Bunch::get_z_period_length)
                .def("get_longitudinal_aperture_length", &Bunch::get_longitudinal_aperture_length)
                .def("get_bucket_barrier_length", &Bunch::get_bucket_barrier_length)
                .def("get_longitudinal_boundary", &Bunch::get_longitudinal_boundary)
                .def("get_bucket_index", &Bunch::get_bucket_index)
                .def("is_periodic", &Bunch::is_z_periodic)
                .def("has_longitudinal_aperture", &Bunch::has_longitudinal_aperture)
                .def("has_bucket_barrier", &Bunch::has_bucket_barrier)
                .def("is_open_longitudinal_boundary", &Bunch::is_open_longitudinal_boundary)
                .def("is_bucket_index_assigned", &Bunch::is_bucket_index_assigned)
                .def("get_comm", &Bunch::get_comm_sptr)
                .def("inject", &Bunch::inject)
                .def("read_file", &Bunch::read_file)
                ;

    enum_<Bunch::State > ("State")
        .value("fixed_z_lab", Bunch::fixed_z_lab)
        .value("fixed_t_lab", Bunch::fixed_t_lab)
        .value("fixed_z_bunch", Bunch::fixed_z_bunch)
        .value("fixed_t_bunch", Bunch::fixed_t_bunch)
        .export_values();

    scope().attr("x") = Bunch::x;
    scope().attr("xp") = Bunch::xp;
    scope().attr("y") = Bunch::y;
    scope().attr("yp") = Bunch::yp;
    scope().attr("z") = Bunch::z;
    scope().attr("zp") = Bunch::zp;
    scope().attr("cdt") = Bunch::cdt;
    scope().attr("dpop") = Bunch::dpop;
    scope().attr("id") = Bunch::id;
}
