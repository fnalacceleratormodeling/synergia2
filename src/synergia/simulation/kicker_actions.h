#ifndef KICKER_ACTIONS_H_
#define KICKER_ACTIONS_H_

#include "synergia/simulation/propagate_actions.h"
#include "synergia/lattice/lattice_element.h"






struct Kick_element
  {   
     Kick_element(std::string const& type, std::string const& name);
     Lattice_element element;
     std::map< int, std::list<int> > map_turn_bunches;    // turns are from 0 to num_turns-1, bunches from 0 to num_bunches-1  
  };
  
typedef   std::list<Kick_element> Kick_elements;

class Kicker_actions : public Propagate_actions
{
  private:
   boost::shared_ptr<Stepper > kstepper_sptr;
   Kick_elements elements_to_kick;
   std::map<std::string, Kick_elements > map_step_to_elements;
   std::list<int> kick_turns;
   void determine_map_step_to_elements();
 public: 
  Kicker_actions(boost::shared_ptr<Stepper > stepper_sptr);
  
  void 
  add_element_to_kick(Kick_element  element);
  
  std::map<std::string, Kick_elements > &
  get_map_step_to_elements() ;
 
  std::list<int> & 
  get_kick_turns() ;
  
 virtual void 
 operator_action(Stepper & stepper, Step & step, Operator * op,
                   Bunch_train & bunch_train, int turn_num, int step_num, int bunch_num);
};


typedef boost::shared_ptr<Kicker_actions > Kicker_actions_sptr;
#endif /* KICKER_ACTIONS_H_ */
