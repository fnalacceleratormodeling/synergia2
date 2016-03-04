#include <boost/python.hpp>
#include "space_charge_3d_open_hockney.h"
#include "space_charge_2d_open_hockney.h"
#include "interpolate_rectangular_zyx.h"
#include "space_charge_2d_bassetti_erskine.h"
#include "space_charge_rectangular.h"
#include "impedance.h"
#include "wake_field.h"
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

void
(Impedance::*apply_bunch)(Bunch &, double, Step &, int, Logger &) = &Impedance::apply;

void
(Impedance::*apply_bunch_train)(Bunch_train &, double, Step &, int,
        Train_diagnosticss const&, Logger &) = &Impedance::apply;

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(set_fftw_helper_member_overloads12,
//        Space_charge_rectangular::set_fftw_helper, 1, 2);

BOOST_PYTHON_MODULE(collective)
{
    class_<Rectangular_grid_domain, Rectangular_grid_domain_sptr>(
        "Rectangular_grid_domain",
        init<std::vector<double> const&, std::vector<double> const&,
             std::vector<int> const&, bool>())
        .def(init<std::vector<double> const&, std::vector<double> const&,
                  std::vector<int> const&, bool>())
        .def(init<std::vector<double> const&, std::vector<int> const&, bool>())
        .def("get_physical_size", &Rectangular_grid_domain::get_physical_size,
             return_value_policy<copy_const_reference>())
        .def("get_physical_offset",
             &Rectangular_grid_domain::get_physical_offset,
             return_value_policy<copy_const_reference>())
        .def("get_grid_shape", &Rectangular_grid_domain::get_grid_shape,
             return_value_policy<copy_const_reference>())
        .def("get_cell_size", &Rectangular_grid_domain::get_cell_size,
             return_value_policy<copy_const_reference>())
        .def("is_periodic", &Rectangular_grid_domain::is_periodic);

    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_3d_open_hockney",
                init<Commxx_divider_sptr, std::vector<int > >())
                .def(init<Commxx_divider_sptr, std::vector<int >, bool >())
                .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool >())
                .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool, double >())
                .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool,
                        double, bool, double >())
                .def(init<std::vector<int > >())
                .def(init<std::vector<int >, bool >())
                .def(init<std::vector<int >, bool, bool >())
                .def(init<std::vector<int >, bool, bool, double >())
                .def(init<std::vector<int >, bool, bool,
                     double, bool >())
                .def(init<std::vector<int >, bool, bool,
                     double, bool >())
                .def(init<std::vector<int >, bool, bool,
                     double, bool, double >())
                .def(init<Commxx_sptr, std::vector<int > >())
                .def(init<Commxx_sptr, std::vector<int >, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool,
                        double, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool,
                        double, bool, double >())
                .def("set_fixed_domain", &Space_charge_3d_open_hockney::set_fixed_domain)
                .def("apply", &Space_charge_3d_open_hockney::apply)
        ;

    class_<Space_charge_2d_open_hockney, Space_charge_2d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_2d_open_hockney",
                init<Commxx_sptr, std::vector<int > >())
        .def(init<Commxx_divider_sptr, std::vector<int > >())
                .def(init<Commxx_sptr, std::vector<int >, bool >())
        .def(init<Commxx_divider_sptr, std::vector<int >, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool >())
        .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double >())
        .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool, double >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double, bool >())
        .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool, double, bool >())
                .def(init<Commxx_sptr, std::vector<int >, bool, bool, double, bool, double >())
        .def(init<Commxx_divider_sptr, std::vector<int >, bool, bool, double, bool, double >())

                .def("apply", &Space_charge_2d_open_hockney::apply)
                .def("set_files", &Space_charge_2d_open_hockney::set_files)
        ;

    scope Space_charge_2d_bassetti_erskine_scope =
        class_<Space_charge_2d_bassetti_erskine, Space_charge_2d_bassetti_erskine_sptr,
        bases<Collective_operator > >("Space_charge_2d_bassetti_erskine",
                init<>())
        .def("apply", &Space_charge_2d_bassetti_erskine::apply)
        .def("set_longitudinal", &Space_charge_2d_bassetti_erskine::set_longitudinal)
        ;
    scope().attr("longitudinal_uniform") = Space_charge_2d_bassetti_erskine::longitudinal_uniform;
    scope().attr("longitudinal_gaussian") = Space_charge_2d_bassetti_erskine::longitudinal_gaussian;

   class_<Space_charge_rectangular, Space_charge_rectangular_sptr,
        bases<Collective_operator > >("Space_charge_rectangular",
              init<Commxx_sptr, std::vector<double >, std::vector<int > , bool >())
           //   .def(init<Commxx_sptr, std::vector<double >, std::vector<int >, bool >())
              .def(init<std::vector<double >, std::vector<int > >())
              //.def("set_fftw_helper", &Space_charge_rectangular::set_fftw_helper,
		//                  set_fftw_helper_member_overloads12())
	      .def("set_fftw_helper", &Space_charge_rectangular::set_fftw_helper)                  
              .def("get_pipe_size", &Space_charge_rectangular::get_pipe_size)
              .def("get_grid_shape", &Space_charge_rectangular::get_grid_shape)
              .def("apply", &Space_charge_rectangular::apply)
        ;

    class_<Wake_field,Wake_field_sptr>("Wake_field",
	        init<std::string const & , std::string const &  >())
         .def("get_wake_file_name", &Wake_field::get_wake_file_name)
         .def("get_wake_type", &Wake_field::get_wake_type)
         .def("get_xw_lead", &Wake_field::get_xw_lead)
         .def("get_xw_trail", &Wake_field::get_xw_trail)
         .def("get_yw_lead", &Wake_field::get_yw_lead)
         .def("get_yw_trail", &Wake_field::get_yw_trail)
         .def("get_z_wake", &Wake_field::get_z_wake)
         .def("multiply_xw_lead", &Wake_field::multiply_xw_lead)
         .def("multiply_xw_trail", &Wake_field::multiply_xw_trail)
         .def("multiply_yw_lead", &Wake_field::multiply_yw_lead)
         .def("multiply_yw_trail", &Wake_field::multiply_yw_trail)
         .def("multiply_z_wake", &Wake_field::multiply_z_wake)
         ;
 
    class_<Impedance,Impedance_sptr,
        bases<Collective_operator > >("Impedance",
                init<std::string const &,std::string const &, int const  &,  double const &, double const &,
		int const, bool, std::vector<int > >())	
	 .def(init<std::string const &,std::string const &, int const  &,  double const &, int const &,
		int const, bool, std::vector<int >  >())	
        .def("get_orbit_length", &Impedance::get_orbit_length)
        .def("get_bunch_spacing", &Impedance::get_bunch_spacing)
        .def("get_z_grid", &Impedance::get_z_grid)
        .def("set_z_grid", &Impedance::set_z_grid)
        .def("get_wake_field", &Impedance::get_wake_field_sptr)
        .def("get_nstored_turns", &Impedance::get_nstored_turns)
        .def("apply", apply_bunch)
        .def("apply", apply_bunch_train)
        ;

	
	
    def("interpolate_rectangular_zyx", &interpolate_rectangular_zyx);
}
