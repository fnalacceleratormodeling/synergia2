#include "lattice_elements_actions.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/stepper.h"

Kick_element::Kick_element()
{
}

Kick_element::Kick_element(std::string const& type, std::string const& name):
element(type,name)
{
}

template<class Archive>
    void
    Kick_element::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(element);
        ar & BOOST_SERIALIZATION_NVP(map_turn_bunches);
    }

template
void
Kick_element::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Kick_element::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Kick_element::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Kick_element::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
        
Kick_element::~Kick_element()
{
}
        
BOOST_CLASS_EXPORT_IMPLEMENT(Kick_element);       
 



const char Lattice_elements_actions::type_name[] = "lattice_elements_actions";

Lattice_elements_actions::Lattice_elements_actions()
{
}  

Lattice_elements_actions::Lattice_elements_actions(Stepper_sptr stepper_sptr):
Propagate_actions(type_name), kstepper_sptr(stepper_sptr), has_inside_operator_actions(true)
{   
}  

//bool const& 
//Lattice_elements_actions::get_has_inside_operator_actions() const
//{
//  return has_inside_operator_actions;
//}  


void 
Lattice_elements_actions::add_element_to_kick(Kick_element element)
{
  // order bunches in map_turn_bunches
    for (std::map< int, std::list<int> >::iterator it=element.map_turn_bunches.begin();
              it!=element.map_turn_bunches.end(); ++it){
        it->second.sort();
        kick_turns.push_back(it->first);
    }
    kick_turns.unique();
    kick_turns.sort();
    elements_to_kick.push_back(element); 
    determine_map_step_to_elements();
}

void
Lattice_elements_actions::print_actions()
{
  for (Kick_elements::const_iterator eit = elements_to_kick.begin(); 
                           eit !=elements_to_kick.end(); ++eit ){
    
        std::cout<<" element: ("<<eit->element.get_type()<<", "<<eit->element.get_name()
        <<") will be changed to: ("<<eit->element.as_string()<<" ), ";
        
            for (std::map< int, std::list<int> >::const_iterator mit=eit->map_turn_bunches.begin();
                        mit!=eit->map_turn_bunches.end(); ++mit){
                  std::cout<<" {at turn "<<mit->first<<", for bunches: [";
                  for(std::list<int> ::const_iterator bit=mit->second.begin();
                          bit!=mit->second.end(); ++bit ){
                          std::cout<<(*bit)<<" ";
                  }
                  std::cout<<"]} ";
            }
        
        
        std::cout<<std::endl;
        
    }
}  
  


 void 
 Lattice_elements_actions::determine_map_step_to_elements()
 {
      if (elements_to_kick.empty()) {      
        this->map_step_to_elements=std::map<std::string, Kick_elements > ();
        this->kick_turns=std::list<int>();
       return;
      }
      
       
      int step_count = 0;
      for (Steps::const_iterator sit = kstepper_sptr->get_steps().begin(); sit
                      != kstepper_sptr->get_steps().end(); ++sit) {
           ++step_count;
            
           for (Operators::const_iterator oit = (*sit)->get_operators().begin();
                           oit != (*sit)->get_operators().end(); ++oit) {
               
               Kick_elements list_elem;//=std::list<Kick_element >();
               if ((*oit)->get_type()=="independent"){              
                   Independent_operator_sptr iop_sptr=boost::dynamic_pointer_cast<Independent_operator >(*oit);
                   for (Lattice_element_slices::const_iterator it = iop_sptr->get_slices().begin();
                                  it != iop_sptr->get_slices().end(); ++it) {
                       
                        for (Kick_elements::const_iterator eit = elements_to_kick.begin(); 
                           eit !=elements_to_kick.end(); ++eit ){                                                                    
                            if (eit->element.get_name()==(*it)->get_lattice_element().get_name()) 
                                          list_elem.push_back(*eit);   
                                            
                        }   
                          
                  }//Lattice_element_slices
                } 
                // making a string key for ther map, step_num + op name        
                std::stringstream pp;
                pp<<step_count<<(*oit)->get_name();
                std::string  key(pp.str());
                if (!(list_elem.empty())) map_step_to_elements[key]=list_elem;
           } //operators                           
     } //steps                     
 }  

std::map<std::string, Kick_elements > &
Lattice_elements_actions::get_map_step_to_elements() 
{
   return map_step_to_elements;
}  

std::list<int> &  
Lattice_elements_actions::get_kick_turns()
{
  return this->kick_turns;
}  


void
Lattice_elements_actions::lattice_elements_action(Stepper & stepper, Step & step, Operator * op, int step_num, int turn_num, 
                   int bunch_num)
{
    try{
       if (&(*kstepper_sptr)!=&stepper)   throw 
             std::runtime_error("kicker stepper_sptr point to a different object than the propagator stepper "); 
     }
     catch(std::exception const& e){    
           std::cout<<e.what()<< std::endl;   
            MPI_Abort(MPI_COMM_WORLD, 123);
     } 
    
  
    if (find(kick_turns.begin(),kick_turns.end(),turn_num)==kick_turns.end()){ //no kick at this turn
      return;
    }      
    else if (op->get_type()=="independent"){
       std::stringstream pp;
       pp<<step_num<<op->get_name();
       std::string key(pp.str());
       std::map<std::string, Kick_elements >::iterator mit=map_step_to_elements.find(key);
       if (mit != map_step_to_elements.end()){
          for (Kick_elements::const_iterator eit=mit->second.begin();
                                                eit!=mit->second.end();++eit){
            
                std::map< int, std::list<int> >::const_iterator ait=eit->map_turn_bunches.find(turn_num);
                if( ait != eit->map_turn_bunches.end()){ 
                    std::list<int>::const_iterator bunch_find=find(ait->second.begin(),ait->second.end(),bunch_num);
                    if (bunch_find !=ait->second.end()){                     
                        for(std::list<Lattice_element_sptr >::const_iterator 
                                  le_it =stepper.get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
                                  le_it != stepper.get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++le_it){ 
                            if ((*le_it)->get_name()==eit->element.get_name() ){
                                std::map<std::string, double > dattr=eit->element.get_double_attributes();                                                                                
                                            //  std::cout<<"turn: "<<turn_num<<"  bunch: "<<bunch_num<<" before update: ";                                             
                                            //  (*le_it)->print();
                                              for (std::map<std::string, double >::const_iterator dait=dattr.begin();
                                                            dait!=dattr.end();++dait){
                                                    (*le_it)->set_double_attribute(dait->first,dait->second);
                                                } 
                                            //    std::cout<<"turn: "<<turn_num<<"  bunch: "<<bunch_num<<" after update: ";                                               
                                             //   (*le_it)->print();  
                            }
                        }
                        stepper.get_lattice_simulator().update();
                    }                    
                }  
          } 
       }              
   }
}

template<class Archive>
    void
    Lattice_elements_actions::serialize(Archive & ar, const unsigned int version)
    {   
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Propagate_actions);
        ar & BOOST_SERIALIZATION_NVP(kstepper_sptr);
        ar & BOOST_SERIALIZATION_NVP(has_inside_operator_actions);
        ar & BOOST_SERIALIZATION_NVP(elements_to_kick);
        ar & BOOST_SERIALIZATION_NVP(map_step_to_elements);
        ar & BOOST_SERIALIZATION_NVP(kick_turns);    
    }
    

   
template

void
Lattice_elements_actions::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lattice_elements_actions::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lattice_elements_actions::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lattice_elements_actions::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);




Lattice_elements_actions::~Lattice_elements_actions()
{
}  
    
BOOST_CLASS_EXPORT_IMPLEMENT(Lattice_elements_actions); 
