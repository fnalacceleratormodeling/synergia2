#include "propagator.h"
#include "synergia/utils/simple_timer.h"

//void
//Propagator::construct()
//{
//    Lattice_element_slices all_slices;
//    for (Steps::const_iterator s_it = stepper.get_steps().begin(); s_it
//            != stepper.get_steps().end(); ++s_it) {
//        for (Operators::const_iterator o_it = (*s_it)->get_operators().begin(); o_it
//                != (*s_it)->get_operators().end(); ++o_it) {
//            all_slices.splice(all_slices.end(), (*o_it)->get_slices());
//        }
//    }
//    chef_lattice.construct_sliced_beamline(all_slices);
//}

Propagator::Propagator(Stepper_sptr stepper_sptr) :
    stepper_sptr(stepper_sptr)
{
}

// Object_to_sptr_hack is a local struct
struct Object_to_sptr_hack
{
    void
    operator()(void const *) const
    {
    }
};

void
Propagator::propagate(Bunch_with_diagnostics & bunch_with_diagnostics, 
             int num_turns, bool verbose)
{

    double t;
    double t_turn, t_turn1;
    int rank = Commxx().get_rank();
    
    Bunch_sptr bunch_sptr=bunch_with_diagnostics.get_bunch_sptr();
    
    
    
   std::ofstream logfile;
    if (rank == 0) logfile.open("log");
    
    for (int turn = 0; turn < num_turns; ++turn) {
        t_turn=MPI_Wtime();
        if (verbose) {
            if (rank == 0) {
                std::cout << "Propagator: turn " << turn + 1 << "/"
                        << num_turns << std::endl;
            }
       }
         bunch_sptr->get_reference_particle().start_repetition();
         t = simple_timer_current();      
        
        for (Multi_diagnostics::iterator dit = bunch_with_diagnostics.get_per_turn_diagnostics().begin(); dit
                != bunch_with_diagnostics.get_per_turn_diagnostics().end(); ++dit) {
            (*dit)->update_and_write();
        }


 
        t = simple_timer_show(t, "diagnostics-turn");
        int step_count = 0;
        int num_steps = stepper_sptr->get_steps().size();             
       for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
               != stepper_sptr->get_steps().end(); ++it) {
           t = simple_timer_current();
          for (Multi_diagnostics::iterator dit = bunch_with_diagnostics.get_per_step_diagnostics().begin(); dit
                  != bunch_with_diagnostics.get_per_step_diagnostics().end(); ++dit) {
             (*dit)->update_and_write();               
          }
            t = simple_timer_show(t, "diagnostics-step");

            ++step_count;
            if (verbose) {
                if (rank == 0) {
                   std::cout << "Propagator:   step " << step_count << "/"
                           << num_steps <<" s= "<<bunch_sptr->get_reference_particle().get_s()<<" trajectory length="<<bunch_sptr->get_reference_particle().get_trajectory_length()<< std::endl;
                }
            }      
             (*it)->apply(*bunch_sptr); 
            // (*it)->apply(*bunch_sptr, bunch_with_diagnostics.get_full2_diagnostics());  
        }
        
       t_turn1= MPI_Wtime();
       if (rank == 0) {
            logfile<<" turn "<<turn + 1<<" : "<< t_turn1-t_turn<< " \n";
            std::cout<<"  turn "<<turn + 1<<" : "<< t_turn1-t_turn<<std::endl;
            logfile.flush();
         }
    }
    t = simple_timer_current();
     for (Multi_diagnostics::iterator dit = bunch_with_diagnostics.get_per_step_diagnostics().begin(); dit
                   != bunch_with_diagnostics.get_per_step_diagnostics().end(); ++dit) {
              (*dit)->update_and_write();               
     }
    t = simple_timer_show(t, "diagnostics-step");
    
         for (Multi_diagnostics::iterator dit = bunch_with_diagnostics.get_per_turn_diagnostics().begin(); dit
                != bunch_with_diagnostics.get_per_turn_diagnostics().end(); ++dit) {
            (*dit)->update_and_write();
     }

    t = simple_timer_show(t, "diagnostics-turn");
     
     if (rank == 0) logfile.close();    


} 




void
Propagator::propagate(Bunch & bunch, int num_turns,
        Diagnostics & per_step_diagnostics, Diagnostics & per_turn_diagnostics,
        bool verbose)
{
    Multi_diagnostics multi_per_step_diagnostics;
    multi_per_step_diagnostics.append(Diagnostics_sptr(&per_step_diagnostics,
            Object_to_sptr_hack()));
    Multi_diagnostics multi_per_turn_diagnostics;
    multi_per_turn_diagnostics.append(Diagnostics_sptr(&per_turn_diagnostics,
            Object_to_sptr_hack()));
    Bunch_sptr bunch_sptr(&bunch,Object_to_sptr_hack());   
         
    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, multi_per_step_diagnostics, multi_per_turn_diagnostics);     
    propagate( bunch_with_diagnostics, num_turns,verbose); 
    
}

void
Propagator::propagate(Bunch & bunch, int num_turns,
        Multi_diagnostics & per_step_diagnostics,
        Multi_diagnostics & per_turn_diagnostics, bool verbose)
{

    Bunch_sptr bunch_sptr(&bunch,Object_to_sptr_hack());
    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, per_step_diagnostics, per_turn_diagnostics);

    propagate( bunch_with_diagnostics, num_turns,verbose); 

}

/*
void
Propagator::propagate(Bunch_train & bunch_train, int num_turns,
        Diagnostics & per_step_diagnostics, Diagnostics & per_turn_diagnostics,
        bool verbose)
{

    Bunch_sptr bunch_sptr(&bunch,Object_to_sptr_hack());
    Diagnostics_sptr per_step_diagnostics_sptr(&per_step_diagnostics,Object_to_sptr_hack());
    Diagnostics_sptr per_turn_diagnostics_sptr(&per_turn_diagnostics,Object_to_sptr_hack());
    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, per_step_diagnostics_sptr, per_turn_diagnostics_sptr);
    
    propagate( bunch_with_diagnostics, num_turns,verbose);    
    Multi_diagnostics multi_per_step_diagnostics;
    multi_per_step_diagnostics.append(Diagnostics_sptr(&per_step_diagnostics,
            Object_to_sptr_hack()));
    Multi_diagnostics multi_per_turn_diagnostics;
    multi_per_turn_diagnostics.append(Diagnostics_sptr(&per_turn_diagnostics,
            Object_to_sptr_hack()));
    propagate(bunch_train, num_turns, multi_per_step_diagnostics,
            multi_per_turn_diagnostics, verbose);
}*/

/*
void
Propagator::propagate(Bunch_train & bunch_train, int num_turns,
        Multi_diagnostics & per_step_diagnostics,
        Multi_diagnostics & per_turn_diagnostics,
       // std::vector<Multi_diagnostics_sptr > & per_step_diagnosticss,
      //  std::vector<Multi_diagnostics_sptr > & per_turn_diagnosticss,
        bool verbose)
{

    int index=0;
    double t;
    double t_turn, t_turn1;
    int rank = Commxx().get_rank();
    
    std::ofstream logfile;
     if (rank == 0) logfile.open("log");
    
    for (int turn = 0; turn < num_turns; ++turn) {
        t_turn=MPI_Wtime();
     //   if (verbose) {
            if (rank == 0) {
                std::cout << "Propagator: turn " << turn + 1 << "/"
                        << num_turns << std::endl;
            }
     //   }
        bunch_train.get_bunch_sptr(index)->get_reference_particle().start_repetition();
        t = simple_timer_current();      
        
        for (Multi_diagnostics::iterator dit = per_turn_diagnostics.begin(); dit
                != per_turn_diagnostics.end(); ++dit) {
            (*dit)->update_and_write();
        }


 
        t = simple_timer_show(t, "diagnostics-turn");
        int step_count = 0;
        int num_steps = stepper_sptr->get_steps().size();             
       for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
               != stepper_sptr->get_steps().end(); ++it) {
           t = simple_timer_current();
           for (Multi_diagnostics::iterator dit = per_step_diagnostics.begin(); dit
                   != per_step_diagnostics.end(); ++dit) {
              (*dit)->update_and_write();               
           }
            t = simple_timer_show(t, "diagnostics-step");

            ++step_count;
            if (verbose) {
//                 if (rank == 0) {
//                     std::cout << "Propagator:   step " << step_count << "/"
//                             << num_steps <<" s= "<<bunch.get_reference_particle().get_s()<<" trajectory length="<<bunch.get_reference_particle().get_trajectory_length()<< std::endl;
//                 }
            }       
           // (*it)->apply(bunch_train.get_bunch_sptr(index)); 
            // (*it)->apply(bunch, per_step_diagnostics);  
        }
        
       t_turn1= MPI_Wtime();
       if (rank == 0) {
        logfile<<" turn "<<turn + 1<<" : "<< t_turn1-t_turn<< " \n";
        std::cout<<"  turn "<<turn + 1<<" : "<< t_turn1-t_turn<<std::endl;
        logfile.flush();
        }
    }
    t = simple_timer_current();
    for (Multi_diagnostics::iterator it = per_step_diagnostics.begin(); it
            != per_step_diagnostics.end(); ++it) {
        (*it)->update_and_write();
    }
    t = simple_timer_show(t, "diagnostics-step");
    for (Multi_diagnostics::iterator it = per_turn_diagnostics.begin(); it
            != per_turn_diagnostics.end(); ++it) {
           (*it)->update_and_write();
    }
    t = simple_timer_show(t, "diagnostics-turn");
     
     if (rank == 0) logfile.close();    


*/















//     if (bunch_train.get_num_bunches() != per_step_diagnosticss.size()) {
//         throw std::runtime_error(
//                 "Propagator: per_step_diagnosticss must have the same length as the bunch_train");
//     }
//     if (bunch_train.get_num_bunches() != per_turn_diagnosticss.size()) {
//         throw std::runtime_error(
//                 "Propagator: per_turn_diagnosticss must have the same length as the bunch_train");
//     }
//     double t;
//     int rank = Commxx().get_rank();
//     for (int turn = 0; turn < num_turns; ++turn) {
//         if (verbose) {
//             if (rank == 0) {
//                 std::cout << "Propagator: turn " << turn + 1 << "/"
//                         << num_turns << std::endl;
//             }
//         }
//         for (int index = 0; index < bunch_train.get_num_bunches(); ++index) {
//             if (bunch_train.is_on_this_rank(index)) {
//                 bunch_train.get_bunch_sptr(index)->get_reference_particle().start_repetition();
//             }
//         }
//         t = simple_timer_current();
//         for (int index = 0; index < bunch_train.get_num_bunches(); ++index) {
//             if (bunch_train.is_on_this_rank(index)) {
//                 for (Multi_diagnostics::iterator dit =
//                         per_turn_diagnosticss[index]->begin(); dit
//                         != per_turn_diagnosticss[index]->end(); ++dit) {
//                     (*dit)->update_and_write();
//                 }
//             }
//         }
//         t = simple_timer_show(t, "diagnostics-turn");
//         int step_count = 0;
//         int num_steps = stepper_sptr->get_steps().size();
//         for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
//                 != stepper_sptr->get_steps().end(); ++it) {
//             t = simple_timer_current();
//             for (int index = 0; index < bunch_train.get_num_bunches(); ++index) {
//                 if (bunch_train.is_on_this_rank(index)) {
//                     for (Multi_diagnostics::iterator dit =
//                             per_step_diagnosticss[index]->begin(); dit
//                             != per_step_diagnosticss[index]->end(); ++dit) {
//                         (*dit)->update_and_write();
//                     }
//                 }
//             }
//             t = simple_timer_show(t, "diagnostics-step");
// 
//             ++step_count;
//             if (verbose) {
//                 if (rank == 0) {
//                     std::cout << "Propagator:   step " << step_count << "/"
//                             << num_steps << std::endl;
//                 }
//             }
//             (*it)->apply_train(bunch_train);
//         }
//     }
//     t = simple_timer_current();
//     for (int index = 0; index < bunch_train.get_num_bunches(); ++index) {
//         if (bunch_train.is_on_this_rank(index)) {
//             for (Multi_diagnostics::iterator dit =
//                     per_step_diagnosticss[index]->begin(); dit
//                     != per_step_diagnosticss[index]->end(); ++dit) {
//                 (*dit)->update_and_write();
//             }
//         }
//     }
//     t = simple_timer_show(t, "diagnostics-step");
//     for (int index = 0; index < bunch_train.get_num_bunches(); ++index) {
//         if (bunch_train.is_on_this_rank(index)) {
//             for (Multi_diagnostics::iterator dit =
//                     per_turn_diagnosticss[index]->begin(); dit
//                     != per_turn_diagnosticss[index]->end(); ++dit) {
//                 (*dit)->update_and_write();
//             }
//         }
//     }
//     t = simple_timer_show(t, "diagnostics-turn");
//}




Propagator::~Propagator()
{

}
