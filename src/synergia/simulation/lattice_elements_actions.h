#ifndef LATTICE_ELEMENTS_ACTIONS_H_
#define LATTICE_ELEMENTS_ACTIONS_H_

#include "synergia/simulation/propagate_actions.h"
#include "synergia/lattice/lattice_element.h"




struct Kick_element
  {  
     Kick_element(); 
     Kick_element(std::string const& type, std::string const& name);
     Lattice_element element;
     std::map< int, std::list<int> > map_turn_bunches;    // turns are from 0 to num_turns-1, bunches from 0 to num_bunches-1
     template<class Archive>
         void
         serialize(Archive & ar, const unsigned int version);
      virtual   
      ~Kick_element();   
  };
  BOOST_CLASS_EXPORT_KEY(Kick_element);
typedef   std::list<Kick_element> Kick_elements;



class Lattice_elements_actions : public Propagate_actions
{
  private:
   bool has_inside_operator_actions;
   boost::shared_ptr<Stepper > kstepper_sptr;
   Kick_elements elements_to_kick;
   std::map<std::string, Kick_elements > map_step_to_elements;
   std::list<int> kick_turns;
   void determine_map_step_to_elements();
 public: 
  Lattice_elements_actions(); 
  Lattice_elements_actions(boost::shared_ptr<Stepper > stepper_sptr);

  bool const& 
  get_has_inside_operator_actions() const;
  
  void 
  add_element_to_kick(Kick_element  element);
  
  std::map<std::string, Kick_elements > &
  get_map_step_to_elements() ;
 
  std::list<int> & 
  get_kick_turns() ;
  
  void
  print_actions();
  
  virtual void 
  lattice_elements_action(Stepper & stepper, Step & step, Operator * op, int step_num, int turn_num, 
                   int bunch_num); 
  template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
  
   static const char type_name[];    
   virtual
   ~Lattice_elements_actions();
    
};
BOOST_CLASS_EXPORT_KEY(Lattice_elements_actions);

typedef boost::shared_ptr<Lattice_elements_actions > Lattice_elements_actions_sptr;
#endif /* LATTICE_ELEMENTS_ACTIONS_H_ */
