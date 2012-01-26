#include "propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/simulation/standard_diagnostics_actions.h"

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
    Propagate_actions empty_propagate_actions;
    propagate(bunch_with_diagnostics, num_turns, empty_propagate_actions,  verbose);
}


void
Propagator::propagate(Bunch_with_diagnostics & bunch_with_diagnostics,
             int num_turns, Propagate_actions & general_actions, bool verbose)
{


try{
        double t, t_total;
        double t_turn, t_turn1;

        t_total = simple_timer_current();

        int rank = Commxx().get_rank();
        Bunch_sptr bunch_sptr=bunch_with_diagnostics.get_bunch_sptr();

        std::ofstream logfile;
        if (rank == 0) logfile.open("log");
        t = simple_timer_current();
        bunch_with_diagnostics.get_diagnostics_actions_sptr()->first_action(*stepper_sptr, *bunch_sptr);
        t = simple_timer_show(t, "diagnostics_first");
        general_actions.first_action(*stepper_sptr, *bunch_sptr);
        t = simple_timer_show(t, "propagate-general_actions");
        for (int turn = 0; turn < num_turns; ++turn) {
            t_turn=MPI_Wtime();
            if (verbose) {
                if (rank == 0) {
                    std::cout << "Propagator: turn " << turn + 1 << "/"
                            << num_turns << std::endl;
                }
            }
            bunch_sptr->get_reference_particle().start_repetition();
            int step_count = 0;
            int num_steps = stepper_sptr->get_steps().size();
            for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                != stepper_sptr->get_steps().end(); ++it) {
                ++step_count;
                if (verbose>0) {
                    if (rank == 0) {
                    std::cout << "Propagator:   step " << step_count << "/"
                            << num_steps <<" s= "<<bunch_sptr->get_reference_particle().get_s()<<" trajectory length="<<bunch_sptr->get_reference_particle().get_trajectory_length()<< std::endl;
                    }
                }
                (*it)->apply(*bunch_sptr);
                /// apply with diagnostics only for testing purposes
                //(*it)->apply(*bunch_sptr, bunch_with_diagnostics.get_per_step_diagnostics());
                t = simple_timer_current();
                bunch_with_diagnostics.get_diagnostics_actions_sptr()->step_end_action(*stepper_sptr, *(*it), *bunch_sptr,
                    turn, step_count);
                t = simple_timer_show(t, "diagnostics-step");
                general_actions.step_end_action(*stepper_sptr, *(*it),*bunch_sptr,
                        turn, step_count);
                t = simple_timer_show(t, "propagate-general_actions-step");
            }

        t_turn1= MPI_Wtime();
        if (rank == 0) {
                logfile<<" turn "<<turn + 1<<" : "<< t_turn1-t_turn<<"   macroparticles="<<bunch_sptr->get_total_num()<< " \n";
                std::cout<<"  turn "<<turn + 1<<" : "<< t_turn1-t_turn<<"   macroparticles="<<bunch_sptr->get_total_num()<<std::endl;
                logfile.flush();
            }
            t = simple_timer_current();

            bunch_with_diagnostics.get_diagnostics_actions_sptr()->turn_end_action(*stepper_sptr, *bunch_sptr, turn);
            t = simple_timer_show(t, "diagnostics-turn");
            general_actions.turn_end_action(*stepper_sptr, *bunch_sptr, turn);
            t = simple_timer_show(t, "propagate-general_actions-turn");

        }
        if (rank == 0) logfile.close();
        simple_timer_show(t_total, "propagate-total");
    }
catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}




void
Propagator::propagate(Bunch & bunch, int num_turns,
        Generalized_diagnostics & per_step_diagnostics, Generalized_diagnostics & per_turn_diagnostics,
        bool verbose)
{

     Bunch_sptr bunch_sptr(&bunch,Object_to_sptr_hack());

     Standard_diagnostics_actions_sptr diagnostics_actions_sptr(new Standard_diagnostics_actions);

     Generalized_diagnostics_sptr per_step_diagnostics_sptr(&per_step_diagnostics,Object_to_sptr_hack());
     diagnostics_actions_sptr->add_per_step(per_step_diagnostics_sptr);

     Generalized_diagnostics_sptr per_turn_diagnostics_sptr(&per_turn_diagnostics,Object_to_sptr_hack());
     diagnostics_actions_sptr->add_per_turn(per_turn_diagnostics_sptr);

     Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr);
     propagate( bunch_with_diagnostics, num_turns,verbose);

}

void
Propagator::propagate(Bunch & bunch, int num_turns,
        Multi_diagnostics & per_step_diagnostics,
        Multi_diagnostics & per_turn_diagnostics, bool verbose)
{

     Bunch_sptr bunch_sptr(&bunch,Object_to_sptr_hack());
     Standard_diagnostics_actions_sptr diagnostics_actions_sptr(new Standard_diagnostics_actions);

     for (Multi_diagnostics::iterator dit = per_step_diagnostics.begin(); dit
                !=  per_step_diagnostics.end(); ++dit) {
            diagnostics_actions_sptr->add_per_step(*dit);
      }
     for (Multi_diagnostics::iterator dit = per_turn_diagnostics.begin(); dit
                !=  per_turn_diagnostics.end(); ++dit) {
            diagnostics_actions_sptr->add_per_turn(*dit);
      }

    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr);
    propagate( bunch_with_diagnostics, num_turns,verbose);


}

void
Propagator::propagate(Bunch_with_diagnostics_train & bunch_diag_train,
             int num_turns,  bool verbose)
{
         Propagate_actions empty_propagate_actions;
         propagate(bunch_diag_train, num_turns, empty_propagate_actions, verbose);
}


void
Propagator::propagate(Bunch_with_diagnostics_train & bunch_diag_train,
             int num_turns,  Propagate_actions & general_actions, bool verbose)
{

try{
        int rank = Commxx().get_rank();
        double t_turn, t_turn1;

        std::ofstream logfile;
        if (rank == 0) logfile.open("log");
        for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
                    if (bunch_diag_train.is_on_this_rank(index)) {
                        Bunch_sptr bunch_sptr=bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
                        bunch_diag_train.get_bunch_diag_sptr(index)->get_diagnostics_actions_sptr()->first_action(*stepper_sptr, *bunch_sptr);
                        general_actions.first_action(*stepper_sptr, *bunch_sptr);
                    }
        }


        for (int turn = 0; turn < num_turns; ++turn) {
                t_turn=MPI_Wtime();
                if (verbose) {
                    if (rank == 0) {
                        std::cout << "Propagator: turn " << turn + 1 << "/"
                                << num_turns << std::endl;
                    }
                }

                for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
                    if (bunch_diag_train.is_on_this_rank(index)) {
                    bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr()->get_reference_particle().start_repetition();
                }
                }


                int step_count = 0;
                int num_steps = stepper_sptr->get_steps().size();
                for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                    != stepper_sptr->get_steps().end(); ++it) {
                    ++step_count;
                    if (verbose) {
                        if (rank == 0) {
                            std::cout << "Propagator:   step " << step_count << "/" << num_steps <<std::endl;
                        }
                    }
                    (*it)->apply(bunch_diag_train);
                    for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
                    if (bunch_diag_train.is_on_this_rank(index)) {
                            Bunch_sptr bunch_sptr=bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
                            bunch_diag_train.get_bunch_diag_sptr(index)->get_diagnostics_actions_sptr()
                                    ->step_end_action(*stepper_sptr, *(*it), *bunch_sptr, turn, step_count);
                            general_actions.step_end_action(*stepper_sptr,*(*it), *bunch_sptr, turn, step_count);
                    }
                    }
                }
                t_turn1= MPI_Wtime();
                if (rank == 0) {
                    logfile<<" turn "<<turn + 1<<" : "<< t_turn1-t_turn<< " \n";
                    std::cout<<"  turn "<<turn + 1<<" : "<< t_turn1-t_turn<<std::endl;
                    logfile.flush();
                }
                for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
                    if (bunch_diag_train.is_on_this_rank(index)) {
                            Bunch_sptr bunch_sptr=bunch_diag_train.get_bunch_diag_sptr(index)->get_bunch_sptr();
                            bunch_diag_train.get_bunch_diag_sptr(index)->get_diagnostics_actions_sptr()
                                    ->turn_end_action(*stepper_sptr, *bunch_sptr, turn);
                            general_actions.turn_end_action(*stepper_sptr,*bunch_sptr, turn);
                    }
                }
        }
        if (rank == 0) logfile.close();
    }
catch (std::exception const& e) {
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}

void
Propagator::propagate(Bunch & bunch, int num_turns,
        Standard_diagnostics_actions & diagnostics_actions, int verbosity)
{
    Propagate_actions empty_propagate_actions;
    propagate(bunch, num_turns, diagnostics_actions, empty_propagate_actions,
            verbosity);
}

 void
 Propagator::propagate(Bunch & bunch, int num_turns,
         Standard_diagnostics_actions & diagnostics_actions,
         Propagate_actions & general_actions, int verbose)
    {

        Bunch_sptr bunch_sptr(&bunch,Object_to_sptr_hack());
        Standard_diagnostics_actions_sptr diagnostics_actions_sptr(&diagnostics_actions,Object_to_sptr_hack());
        Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr, diagnostics_actions_sptr);

        propagate( bunch_with_diagnostics, num_turns,  general_actions, verbose);

 }



Propagator::~Propagator()
{

}
