#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/impedance.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/operator.h"
#include "bunches_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-11;


double const  orbit_length=165.;
double const bunch_spacing=5.;
int const  zgrid=400;


 BOOST_AUTO_TEST_CASE(test_constructor)
 {
 Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid, orbit_length, bunch_spacing,10);
 }
 
BOOST_AUTO_TEST_CASE(test_copy)
 {
 Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid, orbit_length, bunch_spacing,10);
 Impedance imped1(imped);
 imped1.set_z_grid(30);
 
//  std::cout<<" imped zgrid="<<imped.get_z_grid()<<std::endl;
//  std::cout<<" imped wake file="<<imped.get_wake_field_sptr()->get_wake_file_name()<<std::endl;
//  
//  std::cout<<" imped1 zgrid="<<imped1.get_z_grid()<<std::endl;
//  std::cout<<" imped1 wake file="<<imped1.get_wake_field_sptr()->get_wake_file_name()<<std::endl;
 
 Impedance_sptr imped_sptr=Impedance_sptr(new Impedance("test_wake_pp.dat", "XLXTYLYTZpp",
							zgrid, orbit_length, bunch_spacing,100));
 Impedance imped2(*imped_sptr);
 } 


 
 BOOST_FIXTURE_TEST_CASE(test_apply,Bunches_fixture)
 {
       
      Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid, orbit_length, bunch_spacing,3);   
      double step_length=2.3;
      double time_step=10.;
      const int verbosity = 4;
      Step step(step_length);
      Logger logger(0);
      Bunch_train bunch_train(bunches, bunch_spacing);
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),1);      
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),2);
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),3);
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      // 4 times impedance was applied
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),3);
      std::list< std::vector<Bunch_properties> >::const_iterator it=imped.get_stored_vbunches().begin();
      BOOST_CHECK_EQUAL(it->size(),bunches.size()); 
      ++it; ++it;
      BOOST_CHECK_EQUAL(it->size(),bunches.size());
      
 }

BOOST_FIXTURE_TEST_CASE(test_apply_full_machine,Bunches_fixture)
 {
       
      std::vector<int > wn(3,0);
      Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid, orbit_length, bunch_spacing,3,1,wn);   
      double step_length=2.3;
      double time_step=10.;
      const int verbosity = 4;
      Step step(step_length);
      Logger logger(0);
      Bunch_train bunch_train(bunches, bunch_spacing);
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),1);      
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),2);
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),3);
      imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
      // 4 times impedance was applied
      BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),3);
      std::list< std::vector<Bunch_properties> >::const_iterator it=imped.get_stored_vbunches().begin();
      BOOST_CHECK_EQUAL(it->size(),bunches.size()); 
      ++it; ++it;
      BOOST_CHECK_EQUAL(it->size(),bunches.size());      
 }

