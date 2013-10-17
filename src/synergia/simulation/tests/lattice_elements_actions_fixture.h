#ifndef LATTICE_ELEMENTS_ACTIONS_FIXTURE_H_
#define LATTICE_ELEMENTS_ACTIONS_FIXTURE_H_
#include "lattice_fixture.h"
#include "synergia/simulation/lattice_elements_actions.h"
#include "propagator_fixture.h"

namespace kludge
{
#include "synergia/bunch/tests/bunches_fixture.h"
}

struct Lattice_elements_action_fixture
{
  Lattice_elements_action_fixture():
  kick_actions(),
  covariances(boost::extents[6][6]), means(boost::extents[6]), seed(67), distribution(seed, *l.b.comm_sptr)
  {
    
     // kick f at turn1 for bunches 0 and 1   
    Kick_element ef1("quadrupole", "f");   
    ef1.element.set_double_attribute("k1", 0.33);
    
    std::list<int> ef1bunches_for_turn;
    ef1bunches_for_turn.push_back(0);
    ef1bunches_for_turn.push_back(1);
    ef1.map_turn_bunches[1]=ef1bunches_for_turn;
    kick_actions.add_element_to_kick(ef1);
    
    // kick back f at turn1 for bunch 2   
       
    Kick_element ef2("quadrupole", "f"); 
    ef2.element.set_double_attribute("k1", 0.07);
    std::list<int> ef2bunches_for_turn;    
    ef2bunches_for_turn.push_back(2);
    ef2.map_turn_bunches[1]=ef2bunches_for_turn;
    kick_actions.add_element_to_kick(ef2);
    
//  kick back f at turn2 for bunch 0,1,2   
    Kick_element ef3("quadrupole", "f"); 
    ef3.element.set_double_attribute("k1", 0.07);

    std::list<int> ef3bunches_for_turn;
    ef3bunches_for_turn.push_back(0);
    ef3bunches_for_turn.push_back(1);
    ef3bunches_for_turn.push_back(2);
    ef3.map_turn_bunches[2]=ef3bunches_for_turn;
    kick_actions.add_element_to_kick(ef3);
    

    
    Kick_element ed("quadrupole", "d");
    ed.element.set_double_attribute("k1", 0.2);
    ed.element.set_double_attribute("l", 0.7);
    std::list<int> edbunches_for_turn3;
    edbunches_for_turn3.push_back(2);
    edbunches_for_turn3.push_back(0);  
    ed.map_turn_bunches[3]=ef1bunches_for_turn;
    ed.map_turn_bunches[2]=edbunches_for_turn3;
    kick_actions.add_element_to_kick(ed);
    
    
    for (int i = 0; i < 6; ++i) {
            means[i] = 0.; // i * 0.0000072;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1)
                        * 0.00000000001;
            }
            covariances[i][i] *= 10.0; // this makes for a positive-definite matrix
        }

        for (int i = 0; i < 6; ++i) {
            covariances[1][i] *= 0.01;
            covariances[i][1] *= 0.01;
            covariances[3][i] *= 0.01;
            covariances[i][3] *= 0.01;
            covariances[5][i] *= 0.01;
            covariances[i][5] *= 0.01;
        }
    
    
    
  } 
  
  Lattice_fixture l;
  Propagator_fixture p;
  Lattice_elements_actions kick_actions;
  kludge::Bunches_fixture bs;
  MArray2d covariances;
  MArray1d means;
  int seed;
  Random_distribution distribution;
};



#endif /* LATTICE_ELEMENTS_ACTIONS_FIXTURE_H_ */
