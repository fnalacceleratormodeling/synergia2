#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/bunch_train.h"
#include "bunches_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double default_tolerance = 1.0e-14;

BOOST_FIXTURE_TEST_CASE(construct1, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);
}

BOOST_FIXTURE_TEST_CASE(construct2, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);
}

BOOST_FIXTURE_TEST_CASE(construct3, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    const double bunch_separation3 = 99.9;
    separations.push_back(bunch_separation3);
    bool caught_error = false;
    try {
        Bunch_train bunch_train(bunches, separations);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(counts_and_offsets, Bunches_fixture)
{
   const double bunch_separation = 1.7;
   Bunch_train bunch_train(bunches, bunch_separation);
   std::vector<int > counts(bunch_train.get_proc_counts_for_impedance());
   std::vector<int > offsets(bunch_train.get_proc_offsets_for_impedance());
   int num_procs=bunch_train.get_parent_comm_sptr()->get_size();
   int rank=bunch_train.get_parent_comm_sptr()->get_rank();  
   Bunches bunches(bunch_train.get_bunches());
   size_t num_bunches=bunches.size();

// std::cout<<" num_procs="<<num_procs<<"   rank="<<rank<<std::endl;
//     for (int j=0; j < num_bunches; ++j){
//       std::cout<<" bunch["<<j<<"] on this rank "<<bunches[j]->get_comm().has_this_rank()<<"   rank="<<rank<<std::endl;
//       if (bunches[j]->get_comm().has_this_rank())
//                  std::cout<<" bunch["<<j<<"] local_rannk="<<bunches[j]->get_comm().get_rank()<<"   rank="<<rank<<std::endl;
//     }
//  
//     for (int i=0; i < num_procs; ++i){
//      std::cout<<" counts["<<i<<"]="<<counts[i]<<"   rank="<<rank<<std::endl;
//      std::cout<<" ofsets["<<i<<"]="<<offsets[i]<<"   rank="<<rank<<std::endl;     
//    }

  if ((num_bunches==3) && (num_procs==2)) {      
	  BOOST_CHECK_EQUAL(counts[0],1);
	  BOOST_CHECK_EQUAL(offsets[0],0); 
	  BOOST_CHECK_EQUAL(counts[1],2);
	  BOOST_CHECK_EQUAL(offsets[1],1); 
	          	 
  }  
  if ((num_bunches==3) && (num_procs==3)) {       
      BOOST_CHECK_EQUAL(counts[0],1);
      BOOST_CHECK_EQUAL(offsets[0],0); 
      BOOST_CHECK_EQUAL(counts[1],1);
      BOOST_CHECK_EQUAL(offsets[1],1); 
      BOOST_CHECK_EQUAL(counts[2],1);
      BOOST_CHECK_EQUAL(offsets[2],2); 	          	 
  }  
  if ((num_bunches==3) && (num_procs==4)) {
      BOOST_CHECK_EQUAL(counts[0],1);
      BOOST_CHECK_EQUAL(offsets[0],0); 
      BOOST_CHECK_EQUAL(counts[1],1);
      BOOST_CHECK_EQUAL(offsets[1],1); 
      BOOST_CHECK_EQUAL(counts[2],1);
      BOOST_CHECK_EQUAL(offsets[2],2); 
      BOOST_CHECK_EQUAL(counts[3],0);
      BOOST_CHECK_EQUAL(offsets[3],0); 
	  
  }  
  
}

void
compare_bunches(Bunch &bunch1, Bunch &bunch2, double tolerance = default_tolerance,
        bool check_state = true, bool check_ids = true)
{
    BOOST_CHECK_EQUAL(bunch1.get_reference_particle().get_total_energy(),
            bunch2.get_reference_particle().get_total_energy());
    BOOST_CHECK_EQUAL(bunch1.get_particle_charge(),
            bunch2.get_particle_charge());
    BOOST_CHECK_CLOSE(bunch1.get_mass(), bunch2.get_mass(), tolerance);
    BOOST_CHECK_CLOSE(bunch1.get_real_num(), bunch1.get_real_num(), tolerance);
    BOOST_CHECK_EQUAL(bunch1.get_local_num(), bunch2.get_local_num());
    BOOST_CHECK_EQUAL(bunch1.get_total_num(), bunch2.get_total_num());
    BOOST_CHECK_EQUAL(bunch1.get_bucket_index(), bunch2.get_bucket_index());
    if (check_state) {
        BOOST_CHECK_EQUAL(bunch1.get_state(), bunch2.get_state());
    }
    for (int part = 0; part < bunch1.get_local_num(); ++part) {
        // this loop is unrolled in order to give more meaningful error messages
        // i.e., error messages including which component was tested
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][0],
                bunch2.get_local_particles()[part][0], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][1],
                bunch2.get_local_particles()[part][1], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][2],
                bunch2.get_local_particles()[part][2], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][3],
                bunch2.get_local_particles()[part][3], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][4],
                bunch2.get_local_particles()[part][4], tolerance);
        BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][5],
                bunch2.get_local_particles()[part][5], tolerance);
        if (check_ids) {
            BOOST_CHECK_CLOSE(bunch1.get_local_particles()[part][6],
                    bunch2.get_local_particles()[part][6], tolerance);
        }
    }
}


BOOST_FIXTURE_TEST_CASE(serialize_xml, Bunches_fixture)
{
  
   const double bunch_separation = 1.7;
   Bunch_train bunch_train(bunches, bunch_separation);
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   std::string ss("bunch_train");
   std::stringstream pp;
   pp<<rank;
   ss.append(pp.str());
   ss.append(".xml");
   xml_save(bunch_train, ss);
   Bunch_train bunch_loaded;
   xml_load(bunch_loaded, ss);
   
   size_t num_bunches=bunch_train.get_bunches().size();
   size_t num_loaded=bunch_loaded.get_bunches().size();
   BOOST_CHECK_EQUAL(num_bunches, num_loaded);   

   int result;
   MPI_Comm_compare( bunch_loaded.get_parent_comm_sptr()->get(),
 		      bunch_train.get_parent_comm_sptr()->get(), &result);
   BOOST_CHECK(result == MPI_IDENT);  
   
   std::vector<int > counts(bunch_train.get_proc_counts_for_impedance());
   std::vector<int > offsets(bunch_train.get_proc_offsets_for_impedance());
   
   std::vector<int > lcounts(bunch_loaded.get_proc_counts_for_impedance());
   std::vector<int > loffsets(bunch_loaded.get_proc_offsets_for_impedance());
   
   BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(),lcounts.begin(), lcounts.end());
   BOOST_CHECK_EQUAL_COLLECTIONS(offsets.begin(), offsets.end(), loffsets.begin(), loffsets.end());
 
     
     for (int i=0; i<num_bunches; ++i){
       compare_bunches( *bunch_train.get_bunches().at(i), *bunch_loaded.get_bunches().at(i));
     }
    
 
}

