#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/bunch_train.h"
#include "bunches_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

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
