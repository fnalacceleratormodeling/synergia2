#include "propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/simulation/standard_diagnostics_actions.h"
#include "synergia/utils/serialization_files.h"
#include <boost/filesystem.hpp>

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

Propagator::State::State(Bunch_simulator * bunch_simulator_ptr,
        Propagate_actions * propagate_actions_ptr, int num_turns,
        int first_turn, int max_turns, bool verbose) :
    bunch_simulator_ptr(bunch_simulator_ptr),
            propagate_actions_ptr(propagate_actions_ptr), num_turns(num_turns),
            first_turn(first_turn), max_turns(max_turns), verbose(verbose)
{
}

Propagator::Propagator(Stepper_sptr stepper_sptr) :
    stepper_sptr(stepper_sptr), checkpoint_period(10),
            checkpoint_dir(default_checkpoint_dir)
{
}

Propagator::Propagator()
{
}

void
Propagator::set_checkpoint_period(int period)
{
    checkpoint_period = period;
}

int
Propagator::get_checkpoint_period() const
{
    return checkpoint_period;
}

void
Propagator::set_checkpoint_dir(std::string const& directory_name)
{
    checkpoint_dir = directory_name;
}

std::string const&
Propagator::get_checkpoint_dir() const
{
    return checkpoint_dir;
}

// Object_to_sptr_hack is a local struct
struct Object_to_sptr_hack
{
    void
    operator()(void const *) const
    {
    }
};

//void
//Propagator::propagate(Bunch_with_diagnostics & bunch_with_diagnostics,
//        int num_turns, Propagate_actions & general_actions, bool verbose)
void
Propagator::propagate(State & state)
{
    try {
        double t, t_total;
        double t_turn, t_turn1;

        t_total = simple_timer_current();

        int rank = Commxx().get_rank();
        Bunch_sptr bunch_sptr = state.bunch_simulator_ptr->get_bunch_sptr();

        std::ofstream logfile;
        if (rank == 0) logfile.open("log");
        t = simple_timer_current();
        state.bunch_simulator_ptr->get_diagnostics_actions().first_action(
                *stepper_sptr, *bunch_sptr);
        t = simple_timer_show(t, "diagnostics_first");
        state.propagate_actions_ptr->first_action(*stepper_sptr, *bunch_sptr);
        t = simple_timer_show(t, "propagate-general_actions");
        int turns_since_checkpoint = 0;
        int orig_first_turn = state.first_turn;
        for (int turn = state.first_turn; turn < state.num_turns; ++turn) {
            t_turn = MPI_Wtime();
            if (state.verbose) {
                if (rank == 0) {
                    std::cout << "Propagator: turn " << turn + 1 << "/"
                            << state.num_turns << std::endl;
                }
            }
            bunch_sptr->get_reference_particle().start_repetition();
            int step_count = 0;
            int num_steps = stepper_sptr->get_steps().size();
            for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                    != stepper_sptr->get_steps().end(); ++it) {
                ++step_count;
                if (state.verbose > 0) {
                    if (rank == 0) {
                        std::cout << "Propagator:   step " << step_count << "/"
                                << num_steps << " s= "
                                << bunch_sptr->get_reference_particle().get_s()
                                << " trajectory length="
                                << bunch_sptr->get_reference_particle().get_trajectory_length()
                                << std::endl;
                    }
                }
                (*it)->apply(*bunch_sptr);
                /// apply with diagnostics only for testing purposes
                //(*it)->apply(*bunch_sptr, bunch_simulator.get_per_step_diagnostics());
                t = simple_timer_current();
                state.bunch_simulator_ptr->get_diagnostics_actions().step_end_action(
                        *stepper_sptr, *(*it), *bunch_sptr, turn, step_count);
                t = simple_timer_show(t, "diagnostics-step");
                state.propagate_actions_ptr->step_end_action(*stepper_sptr,
                        *(*it), *bunch_sptr, turn, step_count);
                t = simple_timer_show(t, "propagate-general_actions-step");
            }

            t_turn1 = MPI_Wtime();
            if (rank == 0) {
                logfile << " turn " << turn + 1 << " : " << t_turn1 - t_turn
                        << "   macroparticles="
                        << state.bunch_simulator_ptr->get_bunch().get_total_num()
                        << " \n";
                std::cout << "  turn " << turn + 1 << " : " << t_turn1 - t_turn
                        << "   macroparticles="
                        << state.bunch_simulator_ptr->get_bunch().get_total_num()
                        << std::endl;
                logfile.flush();
            }
            t = simple_timer_current();

            state.bunch_simulator_ptr->get_diagnostics_actions().turn_end_action(
                    *stepper_sptr, *bunch_sptr, turn);
            t = simple_timer_show(t, "diagnostics-turn");
            state.propagate_actions_ptr->turn_end_action(*stepper_sptr,
                    *bunch_sptr, turn);
            t = simple_timer_show(t, "propagate-general_actions-turn");
            state.first_turn = turn + 1;
            ++turns_since_checkpoint;
            if ((turns_since_checkpoint == checkpoint_period) && (turn
                    != (state.num_turns - 1))) {
                checkpoint(state);
                turns_since_checkpoint = 0;
            }
            if (((turn - orig_first_turn + 1) == state.max_turns) && (turn
                    != (state.num_turns - 1))) {
                if (rank == 0) {
                    std::cout << "Maximum number of turns reached " << turn
                            << ", " << state.max_turns << std::endl;
                    ;
                }
                if (turns_since_checkpoint > 0) {
                    checkpoint(state);
                }
                break;
            }
        }
        if (rank == 0) logfile.close();
        simple_timer_show(t_total, "propagate-total");
    }
    catch (std::exception const& e) {
        std::cout << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}

const char Propagator::default_checkpoint_dir[] = "checkpoint";

void
Propagator::checkpoint(State & state)
{
    if (state.verbose) {
        if (Commxx().get_rank() == 0) {
            std::cout << "Propagator::checkpoint: ";
        }
    }
    double t0 = MPI_Wtime();
    using namespace boost::filesystem;
    remove_serialization_directory();
    binary_save(*this, get_serialization_path("propagator.bina").c_str());
    binary_save(state, get_serialization_path("state.bina").c_str());
    rename_serialization_directory(checkpoint_dir);
    if (state.verbose) {
        if (Commxx().get_rank() == 0) {
            std::cout << "written to \"" << checkpoint_dir << "\" in "
                    << MPI_Wtime() - t0 << "s\n";
        }
    }
}

Propagator::State
Propagator::get_resume_state(std::string const& checkpoint_directory)
{
    using namespace boost::filesystem;
    State state;
    remove_serialization_directory();
    symlink_serialization_directory(checkpoint_directory);
    binary_load(state,
            get_combined_path(checkpoint_directory, "state.bina").c_str());
    unlink_serialization_directory();
    return state;
}
void
Propagator::resume(std::string const& checkpoint_directory)
{
    State state(get_resume_state(checkpoint_directory));
    propagate(state);
}

void
Propagator::resume(std::string const& checkpoint_directory, int max_turns)
{
    State state(get_resume_state(checkpoint_directory));
    state.max_turns = max_turns;
    propagate(state);
}

void
Propagator::propagate(Bunch_simulator & bunch_simulator, int num_turns,
        int max_turns, bool verbose)
{
    Propagate_actions empty_propagate_actions;
    propagate(bunch_simulator, empty_propagate_actions, num_turns, max_turns,
            verbose);
}

void
Propagator::propagate(Bunch_simulator & bunch_simulator,
        Propagate_actions & general_actions, int num_turns, int max_turns,
        bool verbose)
{
    State state(&bunch_simulator, &general_actions, num_turns, 0, max_turns,
            verbose);
    propagate(state);
}

//void
//Propagator::propagate(Bunch & bunch, int num_turns,
//        Diagnostics & per_step_diagnostics, Diagnostics & per_turn_diagnostics,
//        bool verbose)
//{
//
//    Bunch_sptr bunch_sptr(&bunch, Object_to_sptr_hack());
//
//    Bunch_simulator bunch_simulator(bunch_sptr);
//
//    Diagnostics_sptr per_step_diagnostics_sptr(&per_step_diagnostics,
//            Object_to_sptr_hack());
//    bunch_simulator.get_diagnostics_actions_sptr()->add_per_step(
//            per_step_diagnostics_sptr);
//
//    Diagnostics_sptr per_turn_diagnostics_sptr(&per_turn_diagnostics,
//            Object_to_sptr_hack());
//    bunch_simulator.get_diagnostics_actions_sptr()->add_per_turn(
//            per_turn_diagnostics_sptr);
//
//    propagate(bunch_simulator, num_turns, verbose);
//
//}

//void
//Propagator::propagate(Bunch & bunch, int num_turns,
//        Multi_diagnostics & per_step_diagnostics,
//        Multi_diagnostics & per_turn_diagnostics, bool verbose)
//{
//
//    Bunch_sptr bunch_sptr(&bunch, Object_to_sptr_hack());
//    Standard_diagnostics_actions_sptr diagnostics_actions_sptr(
//            new Standard_diagnostics_actions);
//
//    for (Multi_diagnostics::iterator dit = per_step_diagnostics.begin(); dit
//            != per_step_diagnostics.end(); ++dit) {
//        diagnostics_actions_sptr->add_per_step(*dit);
//    }
//    for (Multi_diagnostics::iterator dit = per_turn_diagnostics.begin(); dit
//            != per_turn_diagnostics.end(); ++dit) {
//        diagnostics_actions_sptr->add_per_turn(*dit);
//    }
//
//    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr,
//            diagnostics_actions_sptr);
//    propagate(bunch_with_diagnostics, num_turns, verbose);
//
//}

void
Propagator::propagate(Bunch_with_diagnostics_train & bunch_diag_train,
        int num_turns, bool verbose)
{
    Propagate_actions empty_propagate_actions;
    propagate(bunch_diag_train, num_turns, empty_propagate_actions, verbose);
}

void
Propagator::propagate(Bunch_with_diagnostics_train & bunch_diag_train,
        int num_turns, Propagate_actions & general_actions, bool verbose)
{

    try {
        int rank = Commxx().get_rank();
        double t_turn, t_turn1;

        std::ofstream logfile;
        if (rank == 0) logfile.open("log");
        for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
            if (bunch_diag_train.is_on_this_rank(index)) {
                Bunch_sptr bunch_sptr = bunch_diag_train.get_bunch_diag_sptr(
                        index)->get_bunch_sptr();
                bunch_diag_train.get_bunch_diag_sptr(index)->get_diagnostics_actions_sptr()->first_action(
                        *stepper_sptr, *bunch_sptr);
                general_actions.first_action(*stepper_sptr, *bunch_sptr);
            }
        }

        for (int turn = 0; turn < num_turns; ++turn) {
            t_turn = MPI_Wtime();
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
                        std::cout << "Propagator:   step " << step_count << "/"
                                << num_steps << std::endl;
                    }
                }
                (*it)->apply(bunch_diag_train);
                for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
                    if (bunch_diag_train.is_on_this_rank(index)) {
                        Bunch_sptr
                                bunch_sptr =
                                        bunch_diag_train.get_bunch_diag_sptr(
                                                index)->get_bunch_sptr();
                        bunch_diag_train.get_bunch_diag_sptr(index)->get_diagnostics_actions_sptr() ->step_end_action(
                                *stepper_sptr, *(*it), *bunch_sptr, turn,
                                step_count);
                        general_actions.step_end_action(*stepper_sptr, *(*it),
                                *bunch_sptr, turn, step_count);
                    }
                }
            }
            t_turn1 = MPI_Wtime();
            if (rank == 0) {
                logfile << " turn " << turn + 1 << " : " << t_turn1 - t_turn
                        << " \n";
                std::cout << "  turn " << turn + 1 << " : " << t_turn1 - t_turn
                        << std::endl;
                logfile.flush();
            }
            for (int index = 0; index < bunch_diag_train.get_num_bunches(); ++index) {
                if (bunch_diag_train.is_on_this_rank(index)) {
                    Bunch_sptr
                            bunch_sptr = bunch_diag_train.get_bunch_diag_sptr(
                                    index)->get_bunch_sptr();
                    bunch_diag_train.get_bunch_diag_sptr(index)->get_diagnostics_actions_sptr() ->turn_end_action(
                            *stepper_sptr, *bunch_sptr, turn);
                    general_actions.turn_end_action(*stepper_sptr, *bunch_sptr,
                            turn);
                }
            }
        }
        if (rank == 0) logfile.close();
    }
    catch (std::exception const& e) {
        std::cout << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}

//void
//Propagator::propagate(Bunch & bunch, int num_turns,
//        Standard_diagnostics_actions & diagnostics_actions, int verbosity)
//{
//    Propagate_actions empty_propagate_actions;
//    propagate(bunch, num_turns, diagnostics_actions, empty_propagate_actions,
//            verbosity);
//}

//void
//Propagator::propagate(Bunch & bunch, int num_turns,
//        Standard_diagnostics_actions & diagnostics_actions,
//        Propagate_actions & general_actions, int verbose)
//{
//
//    Bunch_sptr bunch_sptr(&bunch, Object_to_sptr_hack());
//    Standard_diagnostics_actions_sptr diagnostics_actions_sptr(
//            &diagnostics_actions, Object_to_sptr_hack());
//    Bunch_with_diagnostics bunch_with_diagnostics(bunch_sptr,
//            diagnostics_actions_sptr);
//
//    propagate(bunch_with_diagnostics, num_turns, general_actions, verbose);
//
//}

Propagator::~Propagator()
{

}
