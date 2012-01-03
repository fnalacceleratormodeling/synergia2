#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const int map_order = 1;
const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct_per_elem, Lattice_fixture3)
{
    
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
             "space_charge"));
             
         
     int num_steps_f=1;      
     int num_steps_d=1;
     int num_steps_e=1;
   
   Collective_operators   colective_operators_f;
   colective_operators_f.push_back(space_charge); 
  
   Kicks kicks_f(colective_operators_f, num_steps_f);
   
   Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
             "impedance"));   
   Collective_operators   colective_operators_d;
   colective_operators_d.push_back(impedance);  
   
   Kicks kicks_d(colective_operators_d, num_steps_d);
   
   
   Collective_operators   colective_operators_e; 
   Dummy_collective_operator_sptr dummy(new Dummy_collective_operator(
             "dummy")); 
   colective_operators_e.push_back(dummy); 
  
   Kicks kicks_e(colective_operators_e, num_steps_e);
   
   std::map<std::string, Kicks >  list_choice_map;        
   
   list_choice_map["f"]=kicks_f;
   list_choice_map["d"]=kicks_d;
   list_choice_map["else"]=kicks_e;

  
    bool split_else=false;
     Split_operator_stepper_choice stepper(lattice_simulator, list_choice_map, split_else);
   
    
     BOOST_CHECK_EQUAL(stepper.get_steps().size(),
             lattice_sptr->get_elements().size());
}

BOOST_FIXTURE_TEST_CASE(construct_per_elem1, Lattice_fixture3)
{
    
   
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
             "space_charge"));
             
         
     int num_steps_f=1;      
     int num_steps_d=1;
     int num_steps_e=1;
   
   Collective_operators   colective_operators_f;
   colective_operators_f.push_back(space_charge); 
  
   Kicks kicks_f(colective_operators_f, num_steps_f);
   
   Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
             "impedance"));   
   Collective_operators   colective_operators_d;
   colective_operators_d.push_back(impedance);  
   
   Kicks kicks_d(colective_operators_d, num_steps_d);
   
   
   Collective_operators   colective_operators_e; 
   Dummy_collective_operator_sptr dummy(new Dummy_collective_operator(
             "dummy")); 
   colective_operators_e.push_back(dummy); 
  
   Kicks kicks_e(colective_operators_e, num_steps_e);
   
   std::map<std::string, Kicks >  list_choice_map;        
   
   list_choice_map["f"]=kicks_f;
   list_choice_map["d"]=kicks_d;
   list_choice_map["else"]=kicks_e;

  
    bool split_else=false;
     Split_operator_stepper_choice stepper(lattice_simulator, list_choice_map, split_else);
     double steps_length=0.;
   //  int index=0;
     for (Steps::const_iterator it=stepper.get_steps().begin();it!=stepper.get_steps().end();++it){
         steps_length += (*it)->get_length();
      //   ++index;
       //  (*it)->print(index);
        // std::cout<<" step length="<<(*it)->get_length()<< std::endl;
       //  std::cout<<std::endl;
     }
     // std::cout<<"total step length="<<steps_length<<std::endl; 
      
      double orbit_length=lattice_sptr->get_length();   
         
  
     BOOST_CHECK_CLOSE(orbit_length,  steps_length  , tolerance);    
}         
         
 
BOOST_FIXTURE_TEST_CASE(construct_per_elem2, Lattice_fixture3)
{
    // no "else" keyword
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
             "space_charge"));
             
         
     int num_steps_f=2;      
     int num_steps_d=2;
    
   
   Collective_operators   colective_operators_f;
   colective_operators_f.push_back(space_charge); 
  
   Kicks kicks_f(colective_operators_f, num_steps_f);
   
   Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
             "impedance"));   
   Collective_operators   colective_operators_d;
   colective_operators_d.push_back(impedance);  
   
   Kicks kicks_d(colective_operators_d, num_steps_d);
   

   
   std::map<std::string, Kicks >  list_choice_map;        
   
   list_choice_map["f"]=kicks_f;
   list_choice_map["d"]=kicks_d;
  // list_choice_map["else"]=kicks_e;

  
     int num_steps_e=6;
     bool split_else=false;   
     Split_operator_stepper_choice stepper(num_steps_e, lattice_simulator, list_choice_map,  split_else);
     double steps_length=0.;    
    // int index=0;
     for (Steps::const_iterator it=stepper.get_steps().begin();it!=stepper.get_steps().end();++it){
         steps_length += (*it)->get_length();
       //  ++index;
     //    (*it)->print(index);
     //    std::cout<<" step length="<<(*it)->get_length()<< std::endl;
     //    std::cout<<std::endl;
     }
    //  std::cout<<"total step length="<<steps_length<<std::endl; 
      
      double orbit_length=lattice_sptr->get_length();   
         
  
     BOOST_CHECK_CLOSE(orbit_length,  steps_length  , tolerance);    
    
}        

BOOST_FIXTURE_TEST_CASE(construct_per_elem3, Lattice_fixture3)
{
    // empty list_choice_map
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
//     Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
//              "space_charge"));
             
         
   //  int num_steps_f=2;      
   //  int num_steps_d=2;
    
   
//    Collective_operators   colective_operators_f;
//    colective_operators_f.push_back(space_charge); 
//   
//    Kicks kicks_f(colective_operators_f, num_steps_f);
//    
//    Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
//              "impedance"));   
//    Collective_operators   colective_operators_d;
//    colective_operators_d.push_back(impedance);  
//    
//    Kicks kicks_d(colective_operators_d, num_steps_d);
   

   
   std::map<std::string, Kicks >  list_choice_map;        
   
   //list_choice_map["f"]=kicks_f;
  // list_choice_map["d"]=kicks_d;
  // list_choice_map["else"]=kicks_e;

  
     int num_steps_e=10;
     bool   split_else=false;
     Split_operator_stepper_choice stepper(num_steps_e, lattice_simulator, list_choice_map,  split_else);
     double steps_length=0.;    
  //   int index=0;
     for (Steps::const_iterator it=stepper.get_steps().begin();it!=stepper.get_steps().end();++it){
         steps_length += (*it)->get_length();
     //    ++index;
    //     (*it)->print(index);
    //     std::cout<<" step length="<<(*it)->get_length()<< std::endl;
    //     std::cout<<std::endl;
     }
    //  std::cout<<"total step length="<<steps_length<<std::endl; 
      
      double orbit_length=lattice_sptr->get_length();   
         
  
     BOOST_CHECK_CLOSE(orbit_length,  steps_length  , tolerance);    
    
}        
                  
// BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture3)
// {
//     
//     Lattice_simulator lattice_simulator(lattice_sptr, map_order);
//     Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
//              "space_charge"));
//              
//          
//      int num_steps_f=10;      
//      int num_steps_d=0;
//      int num_steps_e=2;
//    
//    Collective_operators   colective_operators_f;
//    colective_operators_f.push_back(space_charge); 
//   
//    Kicks kicks_f(colective_operators_f, num_steps_f);
//    
//    Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
//              "impedance"));   
//    Collective_operators   colective_operators_d;
//    colective_operators_d.push_back(impedance);  
//    
//    Kicks kicks_d(colective_operators_d, num_steps_d);
//    
//    
//    Collective_operators   colective_operators_e; 
//    Dummy_collective_operator_sptr dummy(new Dummy_collective_operator(
//              "dummy")); 
//    colective_operators_e.push_back(dummy); 
//   
//    Kicks kicks_e(colective_operators_e, num_steps_e);
//    
//    std::map<std::string, Kicks >  list_choice_map;        
//    
//    list_choice_map["f"]=kicks_f;
//    list_choice_map["d"]=kicks_d;
//    list_choice_map["else"]=kicks_e;
//    bool split_else = false;
//      
//          bool caught_error = false;
//      try { 
//         Split_operator_stepper_choice stepper(lattice_simulator, list_choice_map, split_else);
//         
//      }
//      catch (std::runtime_error) {
//          caught_error = true;
//      }
//      BOOST_CHECK(caught_error);
//    
// }         

BOOST_FIXTURE_TEST_CASE(construct_split_else, Lattice_fixture2)
{
    
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
             "space_charge"));
             
         
     int num_steps_f=2;      
     int num_steps_d=1;
     int num_steps_e=6;
   
   Collective_operators   colective_operators_f;
   colective_operators_f.push_back(space_charge); 
  
   Kicks kicks_f(colective_operators_f, num_steps_f);
   
   Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
             "impedance"));   
   Collective_operators   colective_operators_d;
   colective_operators_d.push_back(impedance);  
   
   Kicks kicks_d(colective_operators_d, num_steps_d);
   
   
   Collective_operators   colective_operators_e; 
   Dummy_collective_operator_sptr dummy(new Dummy_collective_operator(
             "dummy")); 
   colective_operators_e.push_back(dummy); 
  
   Kicks kicks_e(colective_operators_e, num_steps_e);
   
   std::map<std::string, Kicks >  list_choice_map;        
   
   list_choice_map["f"]=kicks_f;
   list_choice_map["d"]=kicks_d;
   list_choice_map["else"]=kicks_e;

  
    bool split_else=true;
     Split_operator_stepper_choice stepper(lattice_simulator, list_choice_map, split_else);
         
         double steps_length=0.;    
     // int index=0;
      for (Steps::const_iterator it=stepper.get_steps().begin();it!=stepper.get_steps().end();++it){
          steps_length += (*it)->get_length();
//           ++index;
//           (*it)->print(index);
//           std::cout<<" step length="<<(*it)->get_length()<< std::endl;
//           std::cout<<std::endl;
      } 
      // std::cout<<"total step length="<<steps_length<<std::endl; 

              
       double orbit_length=lattice_sptr->get_length();   
          
   
      BOOST_CHECK_CLOSE(orbit_length,  steps_length  , tolerance);    
    

}

 
BOOST_FIXTURE_TEST_CASE(construct_split_else1, Lattice_fixture4)
{
    
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
             "space_charge"));
             
         
     int num_steps_f=1;      
     int num_steps_d=1;
     int num_steps_e=1;
   
   Collective_operators   colective_operators_f;
   colective_operators_f.push_back(space_charge); 
  
   Kicks kicks_f(colective_operators_f, num_steps_f);
   
   Dummy_collective_operator_sptr impedance(new Dummy_collective_operator(
             "impedance"));   
   Collective_operators   colective_operators_d;
   colective_operators_d.push_back(impedance);  
   
   Kicks kicks_d(colective_operators_d, num_steps_d);
   
   
   Collective_operators   colective_operators_e; 
   Dummy_collective_operator_sptr dummy(new Dummy_collective_operator(
             "dummy")); 
   colective_operators_e.push_back(dummy); 
  
   Kicks kicks_e(colective_operators_e, num_steps_e);
   
   std::map<std::string, Kicks >  list_choice_map;        
   
   list_choice_map["f"]=kicks_f;
   list_choice_map["d"]=kicks_d;
   list_choice_map["else"]=kicks_e;

  
    bool split_else=true;
     Split_operator_stepper_choice stepper(lattice_simulator, list_choice_map, split_else);
         
      double steps_length=0.;    
     // int index=0;
      for (Steps::const_iterator it=stepper.get_steps().begin();it!=stepper.get_steps().end();++it){
          steps_length += (*it)->get_length();
       //   ++index;
        //  (*it)->print(index);
        //  std::cout<<" step length="<<(*it)->get_length()<< std::endl;
        //  std::cout<<std::endl;
      } 
      // std::cout<<"total step length="<<steps_length<<std::endl; 

              
       double orbit_length=lattice_sptr->get_length();   
          
   
      BOOST_CHECK_CLOSE(orbit_length,  steps_length  , tolerance);    
    

}
