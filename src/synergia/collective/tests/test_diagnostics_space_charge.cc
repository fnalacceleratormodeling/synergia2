#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/hdf5_file.h"
#include "gaussian_charge_density.h"
#include "space_charge_bunch_fixtures.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/commxx_divider.h"

 BOOST_GLOBAL_FIXTURE(MPI_fixture);

 const double tolerance = 1.0e-12;


BOOST_FIXTURE_TEST_CASE(update, Ellipsoidal_bunch_fixture)
{
  Diagnostics_space_charge_rectangular diag_spc("spc_rectangular_diagnostics_test.h5");
  Bunch_sptr bunch_sptr((new Bunch(bunch)));
  diag_spc.set_bunch_sptr(bunch_sptr);  
  std::vector<int > grid_shape(3);
  grid_shape[0] = 32;
  grid_shape[1] = 16;
  grid_shape[2] = 18;
  std::vector<double > size(3);
  size[0]=0.1;
  size[1]=0.2;
  size[2]=5.; 
  std::vector<double > offset(3);
  offset[0]=0;
  offset[1]=0;
  offset[2]=0;
  Rectangular_grid grid(size, offset, grid_shape, false);
  int component=1;
  double delta_t=1.e-4;
  double step_beta=12;
  diag_spc.update(bunch, grid, component, delta_t, step_beta);
  diag_spc.write();
  
  
}
BOOST_FIXTURE_TEST_CASE(diagnostics_equally_spread_true, gaussian_3d_bunch_fixture)
{
   Diagnostics_space_charge_rectangular diag_spc("diagnostics_rectangular.h5");
   Diagnostics_space_charge_rectangular_sptr sp_diagnostics_sptr
                  (new Diagnostics_space_charge_rectangular("diagnostics_rectangular.h5"));
   
   
   
  Bunch_sptr bunch_sptr((new Bunch(bunch)));
  bunch_sptr->set_real_num(1.e9);
  MArray1d mean_temp(boost::extents[6]);
  mean_temp=Core_diagnostics::calculate_mean(*bunch_sptr);
  MArray1d std_temp(boost::extents[6]);
  MArray1d mean_bunch(boost::extents[3]);
  mean_bunch[0]= mean_temp[0];
  mean_bunch[1]= mean_temp[2];
  mean_bunch[2]= mean_temp[4];   
  std_temp=Core_diagnostics::calculate_std(*bunch_sptr,mean_temp);
  MArray1d std_bunch(boost::extents[3]);
  std_bunch[0]= std_temp[0];
  std_bunch[1]= std_temp[2];
  std_bunch[2]= std_temp[4];  
 
  
  
  
  sp_diagnostics_sptr->set_bunch_sptr(bunch_sptr);  
   
  double gamma= bunch_sptr->get_reference_particle().get_gamma();
  double beta= bunch_sptr->get_reference_particle().get_beta();
  std::vector<double > size(3);
  size[0]=20*std_bunch[0];
  size[1]=20*std_bunch[1];
  size[2]=10*std_bunch[2];
  std::vector<int > grid_shape(3);  
  grid_shape[0] = 64;
  grid_shape[1] = 128;
  grid_shape[2] = 64;
  
  
  
  
  
  int max_component(3);
  
  double step_length=10.;
  Step step(step_length);
  double betax=20;
  double betay=28;
  double time_step=step_length/beta/pconstants::c;
  step.set_betas(betax, betay);
  
  Space_charge_rectangular space_charge(size, grid_shape);
   
  space_charge.set_fftw_helper(bunch_sptr->get_comm_sptr(), true);
  space_charge.add_diagnostics(sp_diagnostics_sptr);
  space_charge.set_domain(*bunch_sptr);
  Rectangular_grid_sptr rho_sptr(space_charge.get_charge_density(*bunch_sptr));
//  std::vector<double > h(rho_sptr->get_domain_sptr()->get_cell_size());
//   std::ofstream file;
//    file.open("rho.dat");
//    for (int i = 0; i < grid_shape[0]; ++i){
//        file<<"x "<<i*h[0]-0.5*size[0]<<"   "<< rho_sptr->get_grid_points()[i][grid_shape[1]/2-1][grid_shape[2]/2] 
//            <<"   "<< rho_sptr->get_grid_points()[i][grid_shape[1]/2][grid_shape[2]/2] 
//            <<"y "<<i*h[1]-0.5*size[1]<<"   "<<rho_sptr->get_grid_points()[grid_shape[0]/2-1][i][grid_shape[2]/2]
//           <<"   "<<rho_sptr->get_grid_points()[grid_shape[0]/2][i][grid_shape[2]/2] 
//           <<"z "<<i*h[2]-0.5*size[2] <<"   "<<rho_sptr->get_grid_points()[grid_shape[0]/2][grid_shape[1]/2][i]
//           <<std::endl;
//   }  
//      file.close();   
  Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho_sptr,gamma));
  double Enormalization=phi_local->get_normalization();
//    file.open("phi_local.dat");
//    for (int i = 0; i < grid_shape[0]; ++i){
//        file<<"x "<<i*h[0]-0.5*size[0]<<"   "<< phi_local->get_grid_points()[i][grid_shape[1]/2-1][grid_shape[2]/2] 
//            <<"   "<< phi_local->get_grid_points()[i][grid_shape[1]/2][grid_shape[2]/2] 
//            <<"y "<<i*h[1]-0.5*size[1]<<"   "<<phi_local->get_grid_points()[grid_shape[0]/2-1][i][grid_shape[2]/2]
//           <<"   "<<phi_local->get_grid_points()[grid_shape[0]/2][i][grid_shape[2]/2] 
//           <<"z "<<i*h[2]-0.5*size[2] <<"   "<<phi_local->get_grid_points()[grid_shape[0]/2][grid_shape[1]/2][i]
//           <<std::endl;
//   }  
//       file.close();
//   
  
  for (int component = 0; component < max_component; ++component) {
     Rectangular_grid_sptr  En(space_charge.get_En(*phi_local, component));       
//      if (component==0){
//           file.open("Ex.dat");
//           for (int i = 0; i < grid_shape[0]; ++i){
//           file<<"x "<<i*h[0]-0.5*size[0]<<"   "<< En->get_grid_points()[i][grid_shape[1]/2-1][grid_shape[2]/2] 
//               <<"   "<< En->get_grid_points()[i][grid_shape[1]/2][grid_shape[2]/2] 
//               <<"y "<<i*h[1]-0.5*size[1]<<"   "<<En->get_grid_points()[grid_shape[0]/2-1][i][grid_shape[2]/2]
//               <<"   "<<En->get_grid_points()[grid_shape[0]/2][i][grid_shape[2]/2] 
//               <<"z "<<i*h[2]-0.5*size[2] <<"   "<<En->get_grid_points()[grid_shape[0]/2][grid_shape[1]/2][i]
//               <<std::endl;
//           }  
//           file.close();      
//      }
//      else if (component==01){
//           file.open("Ey.dat");
//           for (int i = 0; i < grid_shape[0]; ++i){
//           file<<"x "<<i*h[0]-0.5*size[0]<<"   "<< En->get_grid_points()[i][grid_shape[1]/2-1][grid_shape[2]/2] 
//               <<"   "<< En->get_grid_points()[i][grid_shape[1]/2][grid_shape[2]/2] 
//               <<"y "<<i*h[1]-0.5*size[1]<<"   "<<En->get_grid_points()[grid_shape[0]/2-1][i][grid_shape[2]/2]
//               <<"   "<<En->get_grid_points()[grid_shape[0]/2][i][grid_shape[2]/2] 
//               <<"z "<<i*h[2]-0.5*size[2] <<"   "<<En->get_grid_points()[grid_shape[0]/2][grid_shape[1]/2][i]
//               <<std::endl;
//           }  
//           file.close();             
//      }     
      space_charge.do_diagnostics(*En,component, time_step,step, *bunch_sptr);     
  }
  multi_array_check_equal(sp_diagnostics_sptr->get_mean(),mean_bunch ,  10.e-10);
  multi_array_check_equal(sp_diagnostics_sptr->get_std(),std_bunch ,  10.e-10);
  BOOST_CHECK(floating_point_equal(sp_diagnostics_sptr->get_step_betas()[0], betax, 1e-10));
  BOOST_CHECK(floating_point_equal(sp_diagnostics_sptr->get_step_betas()[1], betay, 1e-10));
  BOOST_CHECK(floating_point_equal(sp_diagnostics_sptr->get_step_length(), step_length, 1e-10));
  
  double Ederiv_x=sp_diagnostics_sptr->get_inc_tune_shift()[0][0]*Enormalization/(sp_diagnostics_sptr->get_factor_tune_Ederiv()*betax);
  double Ederiv_y=sp_diagnostics_sptr->get_inc_tune_shift()[1][0]*Enormalization/(sp_diagnostics_sptr->get_factor_tune_Ederiv()*betay);
  double Q=bunch_sptr->get_real_num()*bunch_sptr->get_particle_charge() * pconstants::e/(std_bunch[2]*std::sqrt(2*mconstants::pi));
    
  double BE_deriv_x=Q/(2*mconstants::pi*pconstants::epsilon0* std_bunch[0]*(std_bunch[0]+std_bunch[1])); //Basseti-Erskin approximation
  double BE_deriv_y=Q/(2*mconstants::pi*pconstants::epsilon0* std_bunch[1]*(std_bunch[0]+std_bunch[1]));//Basseti-Erskin approximation
   
//  std::cout<<" dEx/dx="<<Ederiv_x<<std::endl;
//  std::cout<<" dEy/dy="<<Ederiv_y<<std::endl;       
//  std::cout<<" B dEx/dx="<<BE_deriv_x<<std::endl;
//  std::cout<<" B dEy/dy="<<BE_deriv_y<<std::endl;
   
   BOOST_CHECK(floating_point_equal(Ederiv_x,BE_deriv_x, 5.e-2));
   BOOST_CHECK(floating_point_equal(Ederiv_y,BE_deriv_y, 6.e-2));

}

BOOST_FIXTURE_TEST_CASE(diagnostics_equally_spread_false, gaussian_3d_bunch_fixture)
{  
   Diagnostics_space_charge_rectangular diag_spc("diagnostics_rectangular.h5");
   Diagnostics_space_charge_rectangular_sptr sp_diagnostics_sptr
                  (new Diagnostics_space_charge_rectangular("diagnostics_rectangular.h5"));
   
   
   
  Bunch_sptr bunch_sptr((new Bunch(bunch)));
  bunch_sptr->set_real_num(1.e9);
  MArray1d mean_temp(boost::extents[6]);
  mean_temp=Core_diagnostics::calculate_mean(*bunch_sptr);
  MArray1d std_temp(boost::extents[6]);
  MArray1d mean_bunch(boost::extents[3]);
  mean_bunch[0]= mean_temp[0];
  mean_bunch[1]= mean_temp[2];
  mean_bunch[2]= mean_temp[4];   
  std_temp=Core_diagnostics::calculate_std(*bunch_sptr,mean_temp);
  MArray1d std_bunch(boost::extents[3]);
  std_bunch[0]= std_temp[0];
  std_bunch[1]= std_temp[2];
  std_bunch[2]= std_temp[4];  
 
  
  
  
  sp_diagnostics_sptr->set_bunch_sptr(bunch_sptr);  
   
  double gamma= bunch_sptr->get_reference_particle().get_gamma();
  double beta= bunch_sptr->get_reference_particle().get_beta();
  std::vector<double > size(3);
  size[0]=20*std_bunch[0];
  size[1]=20*std_bunch[1];
  size[2]=10*std_bunch[2];
  std::vector<int > grid_shape(3);  
  grid_shape[0] = 64;
  grid_shape[1] = 128;
  grid_shape[2] = 64;
  
  
  
  int max_component(3);
  
  double step_length=10.;
  Step step(step_length);
  double betax=20;
  double betay=28;
  double time_step=step_length/beta/pconstants::c;
  step.set_betas(betax, betay);
  
  Space_charge_rectangular space_charge(size, grid_shape);
   
  space_charge.set_fftw_helper(bunch_sptr->get_comm_sptr(), false);
  space_charge.add_diagnostics(sp_diagnostics_sptr);
  space_charge.set_domain(*bunch_sptr);
  Rectangular_grid_sptr rho_sptr(space_charge.get_charge_density(*bunch_sptr));
  std::vector<double > h(rho_sptr->get_domain_sptr()->get_cell_size());
  

  Distributed_rectangular_grid_sptr phi_local(space_charge.get_phi_local(*rho_sptr,gamma));
  double Enormalization=phi_local->get_normalization();
 
  std::vector<Rectangular_grid_sptr>  Efield(space_charge.get_Efield(*rho_sptr, *bunch_sptr, max_component,gamma));  
  for (int component = 0; component < max_component; ++component) {
     space_charge.do_diagnostics(*(Efield[component]),component, time_step,step, *bunch_sptr);
   }
  
 
  multi_array_check_equal(sp_diagnostics_sptr->get_mean(),mean_bunch ,  10.e-10);
  multi_array_check_equal(sp_diagnostics_sptr->get_std(),std_bunch ,  10.e-10);
  BOOST_CHECK(floating_point_equal(sp_diagnostics_sptr->get_step_betas()[0], betax, 1e-10));
  BOOST_CHECK(floating_point_equal(sp_diagnostics_sptr->get_step_betas()[1], betay, 1e-10));
  BOOST_CHECK(floating_point_equal(sp_diagnostics_sptr->get_step_length(), step_length, 1e-10));
  
  double Ederiv_x=sp_diagnostics_sptr->get_inc_tune_shift()[0][0]*Enormalization/(sp_diagnostics_sptr->get_factor_tune_Ederiv()*betax);
  double Ederiv_y=sp_diagnostics_sptr->get_inc_tune_shift()[1][0]*Enormalization/(sp_diagnostics_sptr->get_factor_tune_Ederiv()*betay);
  double Q=bunch_sptr->get_real_num()*bunch_sptr->get_particle_charge() * pconstants::e/(std_bunch[2]*std::sqrt(2*mconstants::pi));
    
  double BE_deriv_x=Q/(2*mconstants::pi*pconstants::epsilon0* std_bunch[0]*(std_bunch[0]+std_bunch[1]));//Basseti-Erskin approximation
  double BE_deriv_y=Q/(2*mconstants::pi*pconstants::epsilon0* std_bunch[1]*(std_bunch[0]+std_bunch[1]));//Basseti-Erskin approximation
  
//   std::cout<<" dEx/dx="<<Ederiv_x<<std::endl;
//   std::cout<<" dEy/dy="<<Ederiv_y<<std::endl;  
//   std::cout<<" B dEx/dx="<<BE_deriv_x<<std::endl;
//   std::cout<<" B dEy/dy="<<BE_deriv_y<<std::endl;

  BOOST_CHECK(floating_point_equal(Ederiv_x,BE_deriv_x, 5.e-2));
  BOOST_CHECK(floating_point_equal(Ederiv_y,BE_deriv_y, 6.e-2));

}
