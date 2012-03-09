#include "propagator.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/simulation/diagnostics_actions.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/digits.h"
#include <boost/filesystem.hpp>

const char Propagator::default_checkpoint_dir[] = "checkpoint";
const char Propagator::description_file_name[] = "checkpoint_description.txt";
const char Propagator::propagator_archive_name[] = "propagator.bina";
const char Propagator::state_archive_name[] = "state.bina";
const char Propagator::log_file_name[] = "log";
const char Propagator::stop_file_name[] = "stop";
const char Propagator::alt_stop_file_name[] = "STOPTHISTIMEIMEANIT";

Propagator::State::State(Bunch_simulator * bunch_simulator_ptr,
        Propagate_actions * propagate_actions_ptr, int num_turns,
        int first_turn, int max_turns, int verbosity) :
    bunch_simulator_ptr(bunch_simulator_ptr),
            propagate_actions_ptr(propagate_actions_ptr), num_turns(num_turns),
            first_turn(first_turn), max_turns(max_turns), verbosity(verbosity)
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

void
Propagator::propagate(State & state)
{
    try {
        Logger logger(0, log_file_name);
        double t, t_total;
        double t_turn0, t_turn1;

        t_total = simple_timer_current();
        int rank = Commxx().get_rank();
        Bunch_sptr bunch_sptr = state.bunch_simulator_ptr->get_bunch_sptr();
        t = simple_timer_current();
        if (state.first_turn == 0) {
            state.bunch_simulator_ptr->get_diagnostics_actions().first_action(
                    *stepper_sptr, *bunch_sptr);
            t = simple_timer_show(t, "diagnostics_first");
            state.propagate_actions_ptr->first_action(*stepper_sptr,
                    *bunch_sptr);
            t = simple_timer_show(t, "propagate-general_actions");
        }
        int turns_since_checkpoint = 0;
        int orig_first_turn = state.first_turn;
        bool out_of_particles = false;
        if (state.verbosity > 0) {
            logger << "Propagator: starting turn " << state.first_turn + 1
                    << std::endl;
        }
        for (int turn = state.first_turn; turn < state.num_turns; ++turn) {
            t_turn0 = MPI_Wtime();
            bunch_sptr->get_reference_particle().start_repetition();
            int step_count = 0;
            int num_steps = stepper_sptr->get_steps().size();
            for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it
                    != stepper_sptr->get_steps().end(); ++it) {
                ++step_count;
                double t_step0, t_step1;
                t_step0 = MPI_Wtime();
                (*it)->apply(*bunch_sptr, state.verbosity, logger);
                /// apply with diagnostics only for testing purposes
                //(*it)->apply(*bunch_sptr, bunch_simulator.get_per_step_diagnostics());
                t = simple_timer_current();
                state.bunch_simulator_ptr->get_diagnostics_actions().step_end_action(
                        *stepper_sptr, *(*it), *bunch_sptr, turn, step_count);
                t = simple_timer_show(t, "diagnostics-step");
                state.propagate_actions_ptr->step_end_action(*stepper_sptr,
                        *(*it), *bunch_sptr, turn, step_count);
                t = simple_timer_show(t, "propagate-general_actions-step");
                t_step1 = MPI_Wtime();
                if (state.verbosity > 1) {
                    logger << "Propagator:";
                    logger << "     step " << std::setw(digits(num_steps))
                            << step_count << "/" << num_steps;
                    logger << ", macroparticles = "
                            << state.bunch_simulator_ptr->get_bunch().get_total_num();
                    logger << ", time = " << std::fixed << std::setprecision(3)
                            << t_step1 - t_step0 << "s";
                    ;
                    logger << std::endl;
                }
                if (state.bunch_simulator_ptr->get_bunch().get_total_num() == 0) {
                    if (rank == 0) {
                        std::cout
                                << "Propagator::propagate: No particles left in bunch. Exiting.\n";
                    }
                    out_of_particles = true;
                    break;
                }
            }
            if (out_of_particles) {
                break;
            }
            t = simple_timer_current();

            state.bunch_simulator_ptr->get_diagnostics_actions().turn_end_action(
                    *stepper_sptr, *bunch_sptr, turn);
            t = simple_timer_show(t, "diagnostics-turn");
            state.propagate_actions_ptr->turn_end_action(*stepper_sptr,
                    *bunch_sptr, turn);
            t = simple_timer_show(t, "propagate-general_actions-turn");
            state.first_turn = turn + 1;
            t_turn1 = MPI_Wtime();
            if (state.verbosity > 0) {
                logger << "Propagator:";
                logger << " turn " << std::setw(digits(state.num_turns))
                        << turn + 1 << "/" << state.num_turns;
                logger << ", macroparticles = "
                        << state.bunch_simulator_ptr->get_bunch().get_total_num();
                logger << ", time = " << std::fixed << std::setprecision(2)
                        << t_turn1 - t_turn0 << "s";
                logger << std::endl;
            }
            ++turns_since_checkpoint;
            if ((turns_since_checkpoint == checkpoint_period) && (turn
                    != (state.num_turns - 1))) {
                checkpoint(state, logger, t);
                turns_since_checkpoint = 0;
            }
            if (((turn - orig_first_turn + 1) == state.max_turns) && (turn
                    != (state.num_turns - 1))) {
                logger << "Propagator: maximum number of turns reached\n";
                if (turns_since_checkpoint > 0) {
                    checkpoint(state, logger, t);
                }
                break;
            }
            if (boost::filesystem::exists(stop_file_name)
                    || boost::filesystem::exists(alt_stop_file_name)) {
                logger << "Propagator: stop file detected\n";
                if (turns_since_checkpoint > 0) {
                    checkpoint(state, logger, t);
                }
                break;
            }
        }
        simple_timer_show(t_total, "propagate-total");
    }
    catch (std::exception const& e) {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}

void
Propagator::checkpoint(State & state, Logger & logger, double & t)
{
    if (state.verbosity > 0) {
        logger << "Propagator: checkpoint";
        logger.flush();
    }
    double t0 = MPI_Wtime();
    using namespace boost::filesystem;
    remove_serialization_directory();
    binary_save(*this, get_serialization_path(propagator_archive_name).c_str(),
            true);
    binary_save(state, get_serialization_path(state_archive_name).c_str(), true);
    if (Commxx().get_rank() == 0) {
        std::ofstream description(
                get_serialization_path(description_file_name, false).c_str());
        description << "last_turn = " << state.first_turn << std::endl;
        description << "mpi_size = " << Commxx().get_size() << std::endl;
        description.close();
    }
    rename_serialization_directory(checkpoint_dir);
    double t_checkpoint = MPI_Wtime() - t0;
    if (state.verbosity > 0) {
        logger << " written to \"" << checkpoint_dir << "\"";
        logger << ", time = " << std::fixed << std::setprecision(3)
                << t_checkpoint << "s";
        ;
        logger << std::endl;
    }
    t = simple_timer_show(t, "checkpoint");
}

Propagator::State
Propagator::get_resume_state(std::string const& checkpoint_directory)
{
    using namespace boost::filesystem;
    State state;
    remove_serialization_directory();
    symlink_serialization_directory(checkpoint_directory);
    binary_load(state,
            get_combined_path(checkpoint_directory, state_archive_name).c_str());
    unlink_serialization_directory();
    return state;
}

void
Propagator::resume(std::string const& checkpoint_directory, bool new_max_turns,
        int max_turns, bool new_verbosity, int verbosity)
{
    State state(get_resume_state(checkpoint_directory));
    if (new_max_turns) {
        state.max_turns = max_turns;
    }
    if (new_verbosity) {
        state.verbosity = verbosity;
    }
    propagate(state);
}

void
Propagator::propagate(Bunch_simulator & bunch_simulator, int num_turns,
        int max_turns, int verbosity)
{
    Propagate_actions empty_propagate_actions;
    propagate(bunch_simulator, empty_propagate_actions, num_turns, max_turns,
            verbosity);
}

void
Propagator::propagate(Bunch_simulator & bunch_simulator,
        Propagate_actions & general_actions, int num_turns, int max_turns,
        int verbosity)
{
    State state(&bunch_simulator, &general_actions, num_turns, 0, max_turns,
            verbosity);
    propagate(state);
}

#if 0
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
#endif

Propagator::~Propagator()
{

}
