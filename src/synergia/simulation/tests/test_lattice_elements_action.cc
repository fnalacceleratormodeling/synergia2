#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "propagator_fixture.h"
#include "lattice_elements_actions_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/lattice_elements_actions.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct_lattice_elements_actions, Lattice_elements_action_fixture)
{

    kick_actions.determine_map_step_to_elements(*p.propagator.get_stepper_sptr());
    BOOST_CHECK_EQUAL(kick_actions.get_map_step_to_elements()["1first_half"].front().element.get_name(),"f");
    BOOST_CHECK_EQUAL(kick_actions.get_map_step_to_elements()["4second_half"].front().element.get_name(),"d");
    BOOST_CHECK_EQUAL(kick_actions.get_kick_turns().size(),3);
    BOOST_CHECK_EQUAL(kick_actions.get_kick_turns().front(),1);
    BOOST_CHECK_EQUAL(kick_actions.get_kick_turns().back(),3);
    
    
    //for (std::list<int>::const_iterator it=kick_actions.get_kick_turns().begin();
    //it!=kick_actions.get_kick_turns().end();++it){
    //    std::cout<<" kicked tuns are: "<<*it<<std::endl;
   // }
    
    kick_actions.print_actions();
}  



BOOST_FIXTURE_TEST_CASE(propagate_train_lattice_action, Lattice_elements_action_fixture)
{
   
    const double bunch_spacing = 1.7;
    Bunch_train_sptr bunch_train_sptr(
            new Bunch_train(bs.bunches, bunch_spacing));
    for (Bunches::iterator it = bs.bunches.begin(); it != bs.bunches.end();
            ++it) {
        populate_6d(distribution, **it, means, covariances);
    }
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);


 
    int num_turns = 4;
    p.propagator.propagate(bunch_train_simulator, kick_actions, num_turns);
}




BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_elements_action_fixture)
{

  
    xml_save(kick_actions, "lattice_elements_actions.xml");
    
     Lattice_elements_actions loaded;
     xml_load(loaded, "lattice_elements_actions.xml");
     loaded.print_actions();
}




