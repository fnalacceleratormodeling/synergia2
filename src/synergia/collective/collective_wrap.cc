#include "space_charge_3d_open_hockney.h"
#include "space_charge_2d_open_hockney.h"
#include "interpolate_rectangular_zyx.h"
#include "space_charge_2d_bassetti_erskine.h"
#include "space_charge_rectangular.h"
#include "impedance.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(collective)
{
    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_3d_open_hockney",
                init<Commxx_sptr, std::vector<int > >())
                .def(init<Commxx_sptr, std::vector<int >, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool,
                        double, bool, double >())
                .def("apply", &Space_charge_3d_open_hockney::apply)
                .def("auto_tune_comm", &Space_charge_3d_open_hockney::auto_tune_comm)
        ;

    class_<Space_charge_2d_open_hockney, Space_charge_2d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_2d_open_hockney",
                init<Commxx_sptr, std::vector<int > >())
                .def(init<Commxx_sptr, std::vector<int >, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double, bool, double >())
                .def("apply", &Space_charge_2d_open_hockney::apply)
                .def("set_files", &Space_charge_2d_open_hockney::set_files)
        ;

    class_<Space_charge_2d_bassetti_erskine, Space_charge_2d_bassetti_erskine_sptr,
        bases<Collective_operator > >("Space_charge_2d_bassetti_erskine",
                init<>())
        .def("apply", &Space_charge_2d_bassetti_erskine::apply)
        ;

   class_<Space_charge_rectangular, Space_charge_rectangular_sptr,
        bases<Collective_operator > >("Space_charge_rectangular",
              init<Commxx_sptr, std::vector<double >, std::vector<int > >())
              .def(init<std::vector<double >, std::vector<int > >())
              .def("set_fftw_helper", &Space_charge_rectangular::set_fftw_helper)
              .def("get_pipe_size", &Space_charge_rectangular::get_pipe_size)
              .def("get_grid_shape", &Space_charge_rectangular::get_grid_shape)
              .def("apply", &Space_charge_rectangular::apply)
        ;

    class_<Impedance,Impedance_sptr,
        bases<Collective_operator > >("Impedance",
                init<std::string const &, double const &,  double const &, int const &, std::string const &, int const>())
        .def("get_orbit_length", &Impedance::get_orbit_length)
        .def("get_bunch_spacing", &Impedance::get_bunch_spacing)
        .def("get_z_grid", &Impedance::get_z_grid)
        .def("get_wake_file_name", &Impedance::get_wake_file_name)
        .def("get_nstored_turns", &Impedance::get_nstored_turns)
        .def("apply", &Impedance::apply)
        ;

    def("interpolate_rectangular_zyx", &interpolate_rectangular_zyx);
}
