#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/impedance.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/operator.h"
#include "bunches_fixture.h"
#include "synergia/utils/serialization_files.h"

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

#if 0
BOOST_FIXTURE_TEST_CASE(serialize_, Bunches_fixture)
{
  
  int zgrid1=4;
  std::vector<int > wn(3,0);
  wn[2]=1;
  Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid1, orbit_length, bunch_spacing,3,1,wn);   
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
  
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string ss("imped");
  std::stringstream pp;
  pp<<rank;
  ss.append(pp.str());
  ss.append(".xml");
  xml_save(imped, ss);
  
  Impedance imped_loaded;
  xml_load(imped_loaded, ss);
  
  BOOST_CHECK_EQUAL(imped.is_full_machine() ,  imped_loaded.is_full_machine()  );
  BOOST_CHECK_EQUAL(imped.get_z_grid() ,  imped_loaded.get_z_grid()  );
  BOOST_CHECK_EQUAL(imped.get_orbit_length() ,  imped_loaded.get_orbit_length()  );
  BOOST_CHECK_EQUAL(imped.get_wake_factor() ,  imped_loaded.get_wake_factor()  );
  BOOST_CHECK_EQUAL(imped.get_bunch_spacing() ,  imped_loaded.get_bunch_spacing()  );
  BOOST_CHECK_EQUAL(imped.get_nstored_turns() ,  imped_loaded.get_nstored_turns()  );
  BOOST_CHECK_EQUAL(imped.get_num_buckets() ,  imped_loaded.get_num_buckets() );
  BOOST_CHECK_EQUAL(imped.get_train_wave().size(),imped_loaded.get_train_wave().size());
  for (int i=0;i<imped.get_train_wave().size();++i){
     BOOST_CHECK_EQUAL(imped.get_train_wave().at(i), imped_loaded.get_train_wave().at(i));
  } 
   
 
  std::list< std::vector<Bunch_properties> >::const_iterator it=imped.get_stored_vbunches().begin();
  std::list< std::vector<Bunch_properties> >::const_iterator lit=imped_loaded.get_stored_vbunches().begin();
  for (int st=0; st<imped.get_stored_vbunches().size(); ++st){  
      for (int nb=0; nb<bunches.size();++nb){
	BOOST_CHECK_EQUAL(it->at(nb).x_mean, lit->at(nb).x_mean);
	BOOST_CHECK_EQUAL(it->at(nb).y_mean, lit->at(nb).y_mean);
	BOOST_CHECK_EQUAL(it->at(nb).z_mean, lit->at(nb).z_mean);
	BOOST_CHECK_EQUAL(it->at(nb).realnum, lit->at(nb).realnum);
	BOOST_CHECK_EQUAL(it->at(nb).bucket_index,lit->at(nb).bucket_index);	
      }
      ++it;
      ++lit;
  }
      
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_wake_type(), imped_loaded.get_wake_field_sptr()->get_wake_type());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_wake_file_name(), imped_loaded.get_wake_field_sptr()->get_wake_file_name());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_istart(), imped_loaded.get_wake_field_sptr()->get_istart());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_zstart(), imped_loaded.get_wake_field_sptr()->get_zstart());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_delta_z(), imped_loaded.get_wake_field_sptr()->get_delta_z());
   multi_array_check_equal(imped.get_wake_field_sptr()->get_z_coord(), imped_loaded.get_wake_field_sptr()->get_z_coord(),  tolerance);
   multi_array_check_equal(imped.get_wake_field_sptr()->get_xw_lead(), imped_loaded.get_wake_field_sptr()->get_xw_lead(),  tolerance);
   multi_array_check_equal(imped.get_wake_field_sptr()->get_xw_trail(), imped_loaded.get_wake_field_sptr()->get_xw_trail(),  tolerance);
   multi_array_check_equal(imped.get_wake_field_sptr()->get_yw_lead(), imped_loaded.get_wake_field_sptr()->get_yw_lead(),  tolerance); 
   multi_array_check_equal(imped.get_wake_field_sptr()->get_yw_trail(), imped_loaded.get_wake_field_sptr()->get_yw_trail(),  tolerance); 
   multi_array_check_equal(imped.get_wake_field_sptr()->get_z_wake(), imped_loaded.get_wake_field_sptr()->get_z_wake(),  tolerance);   
   multi_array_check_equal(imped.get_xmom(), imped_loaded.get_xmom(),  tolerance);
   multi_array_check_equal(imped.get_ymom(), imped_loaded.get_ymom(),  tolerance);
   multi_array_check_equal(imped.get_zdensity(), imped_loaded.get_zdensity(),  tolerance);
   multi_array_check_equal(imped.get_xwake_leading(), imped_loaded.get_xwake_leading(),  tolerance);
   multi_array_check_equal(imped.get_xwake_trailing(), imped_loaded.get_xwake_trailing(),  tolerance);
   multi_array_check_equal(imped.get_ywake_leading(), imped_loaded.get_ywake_leading(),  tolerance);
   multi_array_check_equal(imped.get_ywake_trailing(), imped_loaded.get_ywake_trailing(),  tolerance);
   multi_array_check_equal(imped.get_zwake0(), imped_loaded.get_zwake0(),  tolerance);    
  
 
 
}
 
#endif 
 
 BOOST_FIXTURE_TEST_CASE(serialize2_, Bunches_fixture)
{
  
  int zgrid1=4;
  std::vector<int > wn(3,0);
  wn[2]=1;
  Impedance imped("Fwake.dat", "XLXTYLYTZpp",zgrid1, orbit_length, bunch_spacing,3,1,wn);   
  Impedance_sptr imped_clone_sptr(imped.clone()); 
  
  double step_length=2.3;
  double time_step=10.;
  const int verbosity = 1;
  Step step(step_length);
  Logger logger(0);
  Bunch_train bunch_train(bunches, bunch_spacing);
  imped.apply(bunch_train, time_step, step, verbosity, train_diagnosticss, logger);
  BOOST_CHECK_EQUAL(imped.get_stored_vbunches().size(),1);      
  
  imped_clone_sptr->apply(bunch_train, time_step/2, step, verbosity, train_diagnosticss, logger);
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_stored_vbunches().size(),1);      
  
   
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::string ss("imped");
  std::string css("cimped");
  std::stringstream pp;
  pp<<rank;
  ss.append(pp.str());
  ss.append(".xml");
  css.append(pp.str());
  css.append(".xml");
  xml_save(imped, ss);
  xml_save(*imped_clone_sptr, css);
  
  Impedance imped_loaded;
  xml_load(imped_loaded, ss);
 
  Impedance imped_loaded_clone;
  xml_load(imped_loaded_clone, css);
  
  
  BOOST_CHECK_EQUAL(imped.is_full_machine() ,  imped_loaded.is_full_machine()  );
  BOOST_CHECK_EQUAL(imped.get_z_grid() ,  imped_loaded.get_z_grid()  );
  BOOST_CHECK_EQUAL(imped.get_orbit_length() ,  imped_loaded.get_orbit_length()  );
  BOOST_CHECK_EQUAL(imped.get_wake_factor() ,  imped_loaded.get_wake_factor()  );
  BOOST_CHECK_EQUAL(imped.get_bunch_spacing() ,  imped_loaded.get_bunch_spacing()  );
  BOOST_CHECK_EQUAL(imped.get_nstored_turns() ,  imped_loaded.get_nstored_turns()  );
  BOOST_CHECK_EQUAL(imped.get_num_buckets() ,  imped_loaded.get_num_buckets() );
  BOOST_CHECK_EQUAL(imped.get_train_wave().size(),imped_loaded.get_train_wave().size());
  for (int i=0;i<imped.get_train_wave().size();++i){
     BOOST_CHECK_EQUAL(imped.get_train_wave().at(i), imped_loaded.get_train_wave().at(i));
  } 
   
 
  std::list< std::vector<Bunch_properties> >::const_iterator it=imped.get_stored_vbunches().begin();
  std::list< std::vector<Bunch_properties> >::const_iterator lit=imped_loaded.get_stored_vbunches().begin();
  for (int st=0; st<imped.get_stored_vbunches().size(); ++st){  
      for (int nb=0; nb<bunches.size();++nb){
	BOOST_CHECK_EQUAL(it->at(nb).x_mean, lit->at(nb).x_mean);
	BOOST_CHECK_EQUAL(it->at(nb).y_mean, lit->at(nb).y_mean);
	BOOST_CHECK_EQUAL(it->at(nb).z_mean, lit->at(nb).z_mean);
	BOOST_CHECK_EQUAL(it->at(nb).realnum, lit->at(nb).realnum);
	BOOST_CHECK_EQUAL(it->at(nb).bucket_index,lit->at(nb).bucket_index);	
      }
      ++it;
      ++lit;
  }
      
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_wake_type(), imped_loaded.get_wake_field_sptr()->get_wake_type());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_wake_file_name(), imped_loaded.get_wake_field_sptr()->get_wake_file_name());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_istart(), imped_loaded.get_wake_field_sptr()->get_istart());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_zstart(), imped_loaded.get_wake_field_sptr()->get_zstart());
   BOOST_CHECK_EQUAL(imped.get_wake_field_sptr()->get_delta_z(), imped_loaded.get_wake_field_sptr()->get_delta_z());
   multi_array_check_equal(imped.get_wake_field_sptr()->get_z_coord(), imped_loaded.get_wake_field_sptr()->get_z_coord(),  tolerance);
   multi_array_check_equal(imped.get_wake_field_sptr()->get_xw_lead(), imped_loaded.get_wake_field_sptr()->get_xw_lead(),  tolerance);
   multi_array_check_equal(imped.get_wake_field_sptr()->get_xw_trail(), imped_loaded.get_wake_field_sptr()->get_xw_trail(),  tolerance);
   multi_array_check_equal(imped.get_wake_field_sptr()->get_yw_lead(), imped_loaded.get_wake_field_sptr()->get_yw_lead(),  tolerance); 
   multi_array_check_equal(imped.get_wake_field_sptr()->get_yw_trail(), imped_loaded.get_wake_field_sptr()->get_yw_trail(),  tolerance); 
   multi_array_check_equal(imped.get_wake_field_sptr()->get_z_wake(), imped_loaded.get_wake_field_sptr()->get_z_wake(),  tolerance);   
   multi_array_check_equal(imped.get_xmom(), imped_loaded.get_xmom(),  tolerance);
   multi_array_check_equal(imped.get_ymom(), imped_loaded.get_ymom(),  tolerance);
   multi_array_check_equal(imped.get_zdensity(), imped_loaded.get_zdensity(),  tolerance);
   multi_array_check_equal(imped.get_xwake_leading(), imped_loaded.get_xwake_leading(),  tolerance);
   multi_array_check_equal(imped.get_xwake_trailing(), imped_loaded.get_xwake_trailing(),  tolerance);
   multi_array_check_equal(imped.get_ywake_leading(), imped_loaded.get_ywake_leading(),  tolerance);
   multi_array_check_equal(imped.get_ywake_trailing(), imped_loaded.get_ywake_trailing(),  tolerance);
   multi_array_check_equal(imped.get_zwake0(), imped_loaded.get_zwake0(),  tolerance);    

   // now the clones  
  BOOST_CHECK_EQUAL(imped_clone_sptr->is_full_machine() ,  imped_loaded_clone.is_full_machine()  );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_z_grid() ,  imped_loaded_clone.get_z_grid()  );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_orbit_length() ,  imped_loaded_clone.get_orbit_length()  );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_wake_factor() ,  imped_loaded_clone.get_wake_factor()  );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_bunch_spacing() ,  imped_loaded_clone.get_bunch_spacing()  );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_nstored_turns() ,  imped_loaded_clone.get_nstored_turns()  );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_num_buckets() ,  imped_loaded_clone.get_num_buckets() );
  BOOST_CHECK_EQUAL(imped_clone_sptr->get_train_wave().size(),imped_loaded_clone.get_train_wave().size());
  for (int i=0;i<imped_clone_sptr->get_train_wave().size();++i){
     BOOST_CHECK_EQUAL(imped_clone_sptr->get_train_wave().at(i), imped_loaded_clone.get_train_wave().at(i));
  } 
   
 
  it=imped_clone_sptr->get_stored_vbunches().begin();
  lit=imped_loaded_clone.get_stored_vbunches().begin();
  for (int st=0; st<imped_clone_sptr->get_stored_vbunches().size(); ++st){  
      for (int nb=0; nb<bunches.size();++nb){
	BOOST_CHECK_EQUAL(it->at(nb).x_mean, lit->at(nb).x_mean);
	BOOST_CHECK_EQUAL(it->at(nb).y_mean, lit->at(nb).y_mean);
	BOOST_CHECK_EQUAL(it->at(nb).z_mean, lit->at(nb).z_mean);
	BOOST_CHECK_EQUAL(it->at(nb).realnum, lit->at(nb).realnum);
	BOOST_CHECK_EQUAL(it->at(nb).bucket_index,lit->at(nb).bucket_index);	
      }
      ++it;
      ++lit;
  }
      
   BOOST_CHECK_EQUAL(imped_clone_sptr->get_wake_field_sptr()->get_wake_type(), imped_loaded_clone.get_wake_field_sptr()->get_wake_type());
   BOOST_CHECK_EQUAL(imped_clone_sptr->get_wake_field_sptr()->get_wake_file_name(), imped_loaded_clone.get_wake_field_sptr()->get_wake_file_name());
   BOOST_CHECK_EQUAL(imped_clone_sptr->get_wake_field_sptr()->get_istart(), imped_loaded_clone.get_wake_field_sptr()->get_istart());
   BOOST_CHECK_EQUAL(imped_clone_sptr->get_wake_field_sptr()->get_zstart(), imped_loaded_clone.get_wake_field_sptr()->get_zstart());
   BOOST_CHECK_EQUAL(imped_clone_sptr->get_wake_field_sptr()->get_delta_z(), imped_loaded_clone.get_wake_field_sptr()->get_delta_z());
   multi_array_check_equal(imped_clone_sptr->get_wake_field_sptr()->get_z_coord(), imped_loaded_clone.get_wake_field_sptr()->get_z_coord(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_wake_field_sptr()->get_xw_lead(), imped_loaded_clone.get_wake_field_sptr()->get_xw_lead(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_wake_field_sptr()->get_xw_trail(), imped_loaded_clone.get_wake_field_sptr()->get_xw_trail(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_wake_field_sptr()->get_yw_lead(), imped_loaded_clone.get_wake_field_sptr()->get_yw_lead(),  tolerance); 
   multi_array_check_equal(imped_clone_sptr->get_wake_field_sptr()->get_yw_trail(), imped_loaded_clone.get_wake_field_sptr()->get_yw_trail(),  tolerance); 
   multi_array_check_equal(imped_clone_sptr->get_wake_field_sptr()->get_z_wake(), imped_loaded_clone.get_wake_field_sptr()->get_z_wake(),  tolerance);   
   multi_array_check_equal(imped_clone_sptr->get_xmom(), imped_loaded_clone.get_xmom(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_ymom(), imped_loaded_clone.get_ymom(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_zdensity(), imped_loaded_clone.get_zdensity(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_xwake_leading(), imped_loaded_clone.get_xwake_leading(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_xwake_trailing(), imped_loaded_clone.get_xwake_trailing(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_ywake_leading(), imped_loaded_clone.get_ywake_leading(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_ywake_trailing(), imped_loaded_clone.get_ywake_trailing(),  tolerance);
   multi_array_check_equal(imped_clone_sptr->get_zwake0(), imped_loaded_clone.get_zwake0(),  tolerance);
 
} 