
#include "synergia/simulation/propagator.h"

//#include "synergia/utils/simple_timer.h"
#include "synergia/utils/digits.h"

#if 0
// avoid bad interaction between Boost Filesystem and clang
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#endif

#if 0
const std::string Propagator::default_checkpoint_dir = "checkpoint";
const std::string Propagator::description_file_name = "checkpoint_description.txt";
const std::string Propagator::propagator_archive_name = "propagator.bina";
const std::string Propagator::propagator_xml_archive_name = "propagator.xml";
const std::string Propagator::state_archive_name = "state.bina";
const std::string Propagator::state_xml_archive_name = "state.xml";
const std::string Propagator::log_file_name = "log";
const std::string Propagator::stop_file_name = "stop";
const std::string Propagator::alt_stop_file_name = "STOPTHISTIMEIMEANIT";
const int Propagator::default_checkpoint_period = 10;
const int Propagator::default_concurrent_io = 1;
#endif

#if 0
Propagator::State::State(Bunch_simulator * bunch_simulator_ptr,
        Propagate_actions * propagate_actions_ptr, int num_turns,
        int first_turn, int max_turns, int verbosity) :
    bunch_simulator_ptr(bunch_simulator_ptr), bunch_train_simulator_ptr(0),
            propagate_actions_ptr(propagate_actions_ptr), num_turns(num_turns),
            first_turn(first_turn), max_turns(max_turns), verbosity(verbosity)
{
}

Propagator::State::State(Bunch_train_simulator * bunch_train_simulator_ptr,
        Propagate_actions * propagate_actions_ptr, int num_turns,
        int first_turn, int max_turns, int verbosity) :
    bunch_simulator_ptr(0), bunch_train_simulator_ptr(bunch_train_simulator_ptr),
            propagate_actions_ptr(propagate_actions_ptr), num_turns(num_turns),
            first_turn(first_turn), max_turns(max_turns), verbosity(verbosity)
{
}

template<class Archive>
    void
    Propagator::State::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunch_simulator_ptr);
        ar & BOOST_SERIALIZATION_NVP(bunch_train_simulator_ptr);
        ar & BOOST_SERIALIZATION_NVP(propagate_actions_ptr);
        ar & BOOST_SERIALIZATION_NVP(num_turns);
        ar & BOOST_SERIALIZATION_NVP(first_turn);
        ar & BOOST_SERIALIZATION_NVP(max_turns);
        ar & BOOST_SERIALIZATION_NVP(verbosity);
    }

template
void
Propagator::State::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Propagator::State::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Propagator::State::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Propagator::State::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
#endif


#if 0
// omp_threads >= 0 uses the number of threads set by default or by the OMP_NUM_THREADS
// environmental variable
Propagator::Propagator(Stepper_sptr stepper_sptr) :
        stepper_sptr(stepper_sptr), checkpoint_period(
                default_checkpoint_period), checkpoint_dir(
                default_checkpoint_dir), checkpoint_with_xml(false), concurrent_io(
                default_concurrent_io), final_checkpoint(false), omp_threads(-1)
{
}

Propagator::Propagator()
    : omp_threads(-1)
{
}

Stepper_sptr
Propagator::get_stepper_sptr()
{
    return stepper_sptr;
}
#endif

Propagator::Propagator(Lattice const & lattice, Stepper const & stepper)
: lattice(lattice)
, steps()
{
    this->lattice.update();
    steps = stepper.apply(this->lattice);
}

void
Propagator::construct()
{
}

#if 0
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
Propagator::set_checkpoint_with_xml(bool with_xml)
{
    checkpoint_with_xml = with_xml;
}

bool
Propagator::get_checkpoint_with_xml() const
{
    return checkpoint_with_xml;
}

void
Propagator::set_final_checkpoint(bool final_checkpoint)
{
  this->final_checkpoint = final_checkpoint;
}

bool
Propagator::get_final_checkpoint() const
{
  return final_checkpoint;
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

void
Propagator::set_concurrent_io(int max)
{
    concurrent_io = max;
}

int
Propagator::get_concurrent_io() const
{
    return concurrent_io;
}

void
Propagator::set_num_threads(int nt)
{
    if (nt <= 0) {
         return;
    }
    omp_threads = nt;
    omp_set_num_threads(omp_threads);
}

int
Propagator::get_num_threads() const
{
    return omp_threads;
}
#endif

void
Propagator::do_before_start(Bunch_simulator & simulator, Logger & logger)
{
#if 0
    if (state.first_turn == 0) 
    {
        Reference_particle const & lattice_reference_particle = 
                stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle();

        if (state.bunch_simulator_ptr) 
        {
            state.bunch_simulator_ptr->get_diagnostics_actions().first_action(
                    *stepper_sptr, state.bunch_simulator_ptr->get_bunch());

            // set the bunch design_reference_particle from lattice reference particle
            state.bunch_simulator_ptr->get_bunch().set_design_reference_particle(
                    lattice_reference_particle);

        } 
        else 
        {
            size_t num_bunches = state.bunch_train_simulator_ptr->get_bunch_train().get_size();
            for (int i = 0; i < num_bunches; ++i) 
            {
                state.bunch_train_simulator_ptr->get_diagnostics_actionss().at(
                        i)->first_action(*stepper_sptr,
                        *(state.bunch_train_simulator_ptr->get_bunch_train().get_bunches().at(
                                i)));

                // set the bunch design_reference_particle from lattice reference particle
                state.bunch_train_simulator_ptr->get_bunch_train().get_bunches().at(i)
                    ->set_design_reference_particle(lattice_reference_particle);
            }
        }


        t = simple_timer_show(t, "propagate-diagnostics_actions_first");

        if (state.bunch_simulator_ptr) 
        {
            state.propagate_actions_ptr->first_action(*stepper_sptr,
                    state.bunch_simulator_ptr->get_bunch());
        } 
        else 
        {
            state.propagate_actions_ptr->first_action(*stepper_sptr,
                    state.bunch_train_simulator_ptr->get_bunch_train());
        }

        t = simple_timer_show(t, "propagate-general_actions_first");
    }
#endif

    if (simulator.first_turn == 0) 
    {
        simulator.diag_action_step_and_turn(0, 0);
        simulator.prop_action_first(lattice);
        simulator.set_lattice_reference_particle(lattice.get_reference_particle());
    }

    auto updates = lattice.update();
    //if (updates.structure) steps = stepper.apply(lattice);
}

void
Propagator::do_start_repetition(Bunch_simulator & simulator)
{
#if 0
    if (state.bunch_simulator_ptr) {
        state.bunch_simulator_ptr->get_bunch().get_reference_particle().start_repetition();
    } else {
        Bunches & bunches(
                state.bunch_train_simulator_ptr->get_bunch_train().get_bunches());
        for (Bunches::iterator it = bunches.begin(); it != bunches.end();
                ++it) {
            (*it)->get_reference_particle().start_repetition();
        }
    }
#endif

    for(auto & bunch : simulator[0].get_bunches()) 
        bunch.get_reference_particle().start_repetition();

    for(auto & bunch : simulator[1].get_bunches()) 
        bunch.get_reference_particle().start_repetition();
}

void
Propagator::do_step(
        Bunch_simulator & simulator, 
        Step & step, 
        int step_count, 
        int turn_count,
        Logger & logger)
{
    double t_step0 = MPI_Wtime();

    // make sure the lattice is up-to-date
    // e.g., update the chef_lattice after any of the 
    // lattice elements has been updated
    lattice.update();

    // propagate through the step
    step.apply(simulator, logger);

    // t = simple_timer_show(t, "propagate-step_apply");

    // operations associated with bunches
    // e.g., bunch longitudinal operations (periodic, zcut, etc)
    // these operations are not from lattices so are not included 
    // in the steps
    // now I believe this can be included as part of the 
    // prop_action_step_end() method
    // simulator.bunch_operation_step_end();
    // t = simple_timer_show(t, "propagate-bunch_operations_step");

    // general diagnostics
    simulator.diag_action_step_and_turn(turn_count, step_count);
    // t = simple_timer_show(t, "propagate-diagnostics_actions_step");

    // propagate action
    simulator.prop_action_step_end(lattice, turn_count, step_count);
    // t = simple_timer_show(t, "propagate-general_actions-step");

    double t_step1 = MPI_Wtime();

    //if (logger.verbosity > 1) 
    {
        int p = std::cout.precision();
        int num_steps = steps.size();

        logger << "Propagator:     step " << std::setw(digits(num_steps)) 
               << step_count << "/" << num_steps;

        logger << ", s_n=" << std::fixed << std::setprecision(4)
               << simulator[0][0].get_reference_particle().get_s_n();

        logger << ", time = " << std::fixed << std::setprecision(3) 
               << t_step1 - t_step0 << "s" << std::endl;

        std::cout.precision(p);
    }

#if 0
    double t_step0 = MPI_Wtime();

    if (state.bunch_simulator_ptr) 
    {
        Bunch & bunch(state.bunch_simulator_ptr->get_bunch());

        Diagnostics_actions & diagnostics_actions(
                state.bunch_simulator_ptr->get_diagnostics_actions());

        step.apply(bunch, state.verbosity,
                diagnostics_actions.get_per_operator_diagnosticss(),
                diagnostics_actions.get_per_operation_diagnosticss(), *stepper_sptr, logger);       

        t = simple_timer_show(t, "propagate-step_apply");

        diagnostics_actions.step_end_action(
                *stepper_sptr, step, bunch, turn, step_count);

        t = simple_timer_show(t, "propagate-diagnostics_actions_step");

        state.propagate_actions_ptr->step_end_action(
                *stepper_sptr, step, bunch, turn, step_count);
    } 
    else 
    {
        Bunch_train & bunch_train(
                state.bunch_train_simulator_ptr->get_bunch_train());
        Diagnostics_actionss & diagnostics_actionss(
                state.bunch_train_simulator_ptr->get_diagnostics_actionss());
        size_t num_bunches =
                state.bunch_train_simulator_ptr->get_bunch_train().get_size();
        Train_diagnosticss per_operator_train_diagnosticss(num_bunches),
                per_operation_train_diagnosticss(num_bunches);
        for (int i = 0; i < num_bunches; ++i) {
            per_operator_train_diagnosticss.at(i) =
                    diagnostics_actionss.at(i)->get_per_operator_diagnosticss();
            per_operation_train_diagnosticss.at(i) =
                    diagnostics_actionss.at(i)->get_per_operation_diagnosticss();
        }
      
        
          step.apply(bunch_train, state.verbosity, per_operator_train_diagnosticss,
                    per_operation_train_diagnosticss, 
                    state.propagate_actions_ptr, *stepper_sptr, step_count, turn,
                    logger);  
                
                
        t = simple_timer_show(t, "propagate-step_apply");
        for (int i = 0; i < num_bunches; ++i) {
            diagnostics_actionss.at(i)->step_end_action(*stepper_sptr, step,
                    *(state.bunch_train_simulator_ptr->get_bunch_train().get_bunches().at(
                            i)), turn, step_count);
        }
        t = simple_timer_show(t, "propagate-diagnostics_actions_step");
        state.propagate_actions_ptr->step_end_action(*stepper_sptr, step,
                bunch_train, turn, step_count);
    }

    t = simple_timer_show(t, "propagate-general_actions-step");
    double t_step1 = MPI_Wtime();
    if (state.verbosity > 1) {
        int p = cout.precision();
        logger << "Propagator:";
        logger << "     step " << std::setw(digits(num_steps)) << step_count
                << "/" << num_steps;
        if (state.bunch_train_simulator_ptr) {
            logger << ", s_n=" << std::fixed << std::setprecision(4)
                    << state.bunch_train_simulator_ptr->get_bunch_train().get_bunches()[0]->get_reference_particle().get_s_n();

        }
        if (state.bunch_simulator_ptr) {
            logger << ", s_n=" << std::fixed << std::setprecision(4)
                    << state.bunch_simulator_ptr->get_bunch().get_reference_particle().get_s_n();
            logger << ", macroparticles = "
                    << state.bunch_simulator_ptr->get_bunch().get_total_num();
        }
        logger << ", time = " << std::fixed << std::setprecision(3)
                << t_step1 - t_step0 << "s";
        logger << std::endl;
        cout.precision(p);
    }
#endif
}

bool
Propagator::check_out_of_particles(Bunch_simulator const & simulator, Logger & logger) 
{
    return false;
#if 0
    // n.b.: We only check out_of_particles for single-bunch propagation.
    // Checking all bunches in a multi-bunch simulation would require a
    // global communication. The costs exceed the potential benefits.
	bool retval = false;
    if (state.bunch_simulator_ptr) {
        if (state.bunch_simulator_ptr->get_bunch().get_total_num() == 0) {
            logger
                    << "Propagator::propagate: No particles left in bunch. Exiting.\n";
            retval = true;
        }
    }
	return retval;
#endif
}

void
Propagator::do_turn_end(
        Bunch_simulator & simulator, 
        int turn_count, 
		Logger & logger) 
{
    //t = simple_timer_current();

    // diagnostic actions
    simulator.diag_action_step_and_turn(turn_count, Propagator::FINAL_STEP);
    simulator.diag_action_particle_loss_write();

    // propagate actions
    simulator.prop_action_turn_end(lattice, turn_count);

    // update lattice in case it has been chaged in the propagate action
    auto updates = lattice.update();
    //if (updates.structure()) steps = stepper.apply(lattice);

    // increment the turn number
    simulator.first_turn = turn_count + 1;

    //double t_turn1 = MPI_Wtime();

#if 0
    if (state.verbosity > 0) 
    {
        int p = cout.precision();
        logger << "Propagator:";
        logger << " turn " << std::setw(digits(state.num_turns)) << turn + 1
                << "/" << state.num_turns;
        if (state.bunch_simulator_ptr) 
        {
            logger << ", macroparticles = "
                    << state.bunch_simulator_ptr->get_bunch().get_total_num();
        } 
        else 
        {
            Bunch_train & bunch_train(
                    state.bunch_train_simulator_ptr->get_bunch_train());
            size_t num_bunches = bunch_train.get_size();
            bunch_train.update_bunch_total_num();
            Bunches & bunches(bunch_train.get_bunches());
            logger << ", macroparticles = (";
            for (std::vector<Bunch_sptr >::const_iterator bit = bunches.begin(); bit!=bunches.end(); ++bit) {
                if (bit != bunches.begin()) {
                    logger << ", " ;
                }
                logger << (*bit)->get_total_num();
            }
             logger << ")  " ;
        }
        logger << ", time = " << std::fixed << std::setprecision(4)
                << t_turn1 - t_turn0 << "s";
        logger << std::endl;
        cout.precision(p);
    }
#endif

#if 0
    t = simple_timer_current();

    if (state.bunch_simulator_ptr) 
    {
        state.bunch_simulator_ptr->get_diagnostics_actions().turn_end_action(
                *stepper_sptr, state.bunch_simulator_ptr->get_bunch(), turn);
    } 
    else 
    {
        size_t num_bunches =
                state.bunch_train_simulator_ptr->get_bunch_train().get_size();
        for (int i = 0; i < num_bunches; ++i) {
            state.bunch_train_simulator_ptr->get_diagnostics_actionss().at(i)->turn_end_action(
                    *stepper_sptr,
                    *(state.bunch_train_simulator_ptr->get_bunch_train().get_bunches().at(
                            i)), turn);
        }
    }

    t = simple_timer_show(t, "propagate-diagnostics_turn");

    if (state.bunch_simulator_ptr) 
    {
        state.propagate_actions_ptr->turn_end_action(*stepper_sptr,
                state.bunch_simulator_ptr->get_bunch(), turn);
    } 
    else 
    {
        state.propagate_actions_ptr->turn_end_action(*stepper_sptr,
                state.bunch_train_simulator_ptr->get_bunch_train(), turn);
    }

    t = simple_timer_show(t, "propagate-general_actions_turn");

    if (stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_have_loss_diagnostics())
    {
         Diagnostics_losses diagnostics_list =
             stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_loss_diagnostics_list();

         for (Diagnostics_losses::const_iterator d_it = diagnostics_list.begin();
             d_it != diagnostics_list.end(); ++d_it)
         { 
             (*d_it)->write();              
         }
    }
        
    
    
    state.first_turn = turn + 1;
    double t_turn1 = MPI_Wtime();
    if (state.verbosity > 0) {
        int p = cout.precision();
        logger << "Propagator:";
        logger << " turn " << std::setw(digits(state.num_turns)) << turn + 1
                << "/" << state.num_turns;
        if (state.bunch_simulator_ptr) {
            logger << ", macroparticles = "
                    << state.bunch_simulator_ptr->get_bunch().get_total_num();
        } else {
            Bunch_train & bunch_train(
                    state.bunch_train_simulator_ptr->get_bunch_train());
            size_t num_bunches = bunch_train.get_size();
            bunch_train.update_bunch_total_num();
            Bunches & bunches(bunch_train.get_bunches());
            logger << ", macroparticles = (";
            for (auto const & bunch : bunches)
            {
            }
            for (std::vector<Bunch_sptr >::const_iterator bit = bunches.begin(); bit!=bunches.end(); ++bit) {
                if (bit != bunches.begin()) {
                    logger << ", " ;
                }
                logger << (*bit)->get_total_num();
            }
             logger << ")  " ;
        }
        logger << ", time = " << std::fixed << std::setprecision(4)
                << t_turn1 - t_turn0 << "s";
        logger << std::endl;
        cout.precision(p);
    }
#endif
}

void
Propagator::propagate(Bunch_simulator & sim)
{
    try 
    {
        Logger logger(0, "log");

        //double t_total = simple_timer_current();
        //double t = simple_timer_current();

        do_before_start(sim, logger);

        int turns_since_checkpoint = 0;
        int orig_first_turn = sim.first_turn;
        bool out_of_particles = false;

        //if (sim.verbosity > 0)
        {
            logger << "Propagator: starting turn " 
                   << sim.first_turn + 1 << std::endl;
        }

        for (int turn_count = sim.first_turn; turn_count < sim.num_turns; ++turn_count) 
        {
            double t_turn0 = MPI_Wtime();

            do_start_repetition(sim);

            int step_count = 0;
            int num_steps = steps.size();

            for (auto & step : steps)
            {
                ++step_count;
            	do_step(sim, step, step_count, turn_count, logger);

                out_of_particles = check_out_of_particles(sim, logger);
                if (out_of_particles) break;
            }

            if (out_of_particles) break;

            ++turns_since_checkpoint;
            do_turn_end(sim, turn_count, logger);

#if 0
            if ((turns_since_checkpoint == checkpoint_period) 
                    || ((turn_count == (sim.num_turns - 1)) && final_checkpoint)) 
            {
                t = simple_timer_current();
                checkpoint(sim, logger, t);
                t = simple_timer_show(t, "propagate-checkpoint_period");
                turns_since_checkpoint = 0;
            }
#endif

            if (((turn_count - orig_first_turn + 1) == sim.max_turns) 
                    && (turn_count != (sim.num_turns - 1))) 
            {
                logger << "Propagator: maximum number of turns reached\n";

#if 0
                if (turns_since_checkpoint > 0) 
                {
                    t = simple_timer_current();
                    checkpoint(sim, logger, t);
                    t = simple_timer_show(t, "propagate-checkpoint_max");
                }
#endif

                break;
            }

#if 0
            if (boost::filesystem::exists(stop_file_name)
                    || boost::filesystem::exists(alt_stop_file_name)) 
            {
                logger << "Propagator: stop file detected\n";

                if (turns_since_checkpoint > 0) 
                {
                    t = simple_timer_current();
                    checkpoint(sim, logger, t);
                    t = simple_timer_show(t, "propagate-checkpoint_stop");
                }

                break;
            }
#endif

        }

        if (out_of_particles) logger << "Propagator: no particles left\n";
        // simple_timer_show(t_total, "propagate-total");
    }
    catch (std::exception const& e) 
    {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }
}

#if 0
void
Propagator::checkpoint(Bunch_simulator & sim, Logger & logger, double & t)
{
    if (state.verbosity > 0) {
        logger << "Propagator: checkpoint";
        logger.flush();
    }
    double t0 = MPI_Wtime();
    using namespace boost::filesystem;
    remove_serialization_directory();
    Commxx commxx_world;
    const int verbosity_threshold = 2;
    Logger iocclog("iocycle_checkpoint", state.verbosity > verbosity_threshold);
    int max_writers;
    if (concurrent_io == 0) {
        max_writers = commxx_world.get_size();
    } else {
        max_writers = concurrent_io;
    }
    int num_cycles = (commxx_world.get_size() + max_writers - 1) / max_writers;
    for (int cycle = 0; cycle < num_cycles; ++cycle) {
        iocclog << "start cycle " << cycle << std::endl;
        int cycle_min = cycle * max_writers;
        int cycle_max = (cycle + 1) * max_writers;
        if ((commxx_world.get_rank() >= cycle_min)
                && (commxx_world.get_rank() < cycle_max)) {
            iocclog << "start write" << std::endl;
            binary_save(*this,
                    get_serialization_path(propagator_archive_name).c_str(),
                    true);
            binary_save(state,
                    get_serialization_path(state_archive_name).c_str(), true);
            if (checkpoint_with_xml) {
                xml_save(*this,
                        get_serialization_path(propagator_xml_archive_name).c_str(),
                        true);
                xml_save(state,
                        get_serialization_path(state_xml_archive_name).c_str(),
                        true);
            }
            iocclog << "end write" << std::endl;
        }
        MPI_Barrier(commxx_world.get());
    }
    if (commxx_world.get_rank() == 0) {
        std::ofstream description(
                get_serialization_path(description_file_name, false).c_str());
        description << "last_turn=" << state.first_turn << std::endl;
        description << "mpi_size=" << Commxx().get_size() << std::endl;
        description.close();
    }
    rename_serialization_directory(checkpoint_dir);
    double t_checkpoint = MPI_Wtime() - t0;
    if (state.verbosity > 0) {
        int p = cout.precision();
        logger << " written to \"" << checkpoint_dir << "\"";
        logger << ", time = " << std::fixed << std::setprecision(3)
                << t_checkpoint << "s";
        ;
        logger << std::endl;
        cout.precision(p);
    }
    t = simple_timer_show(t, "checkpoint");
}
#endif

#if 0
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
#endif

#if 0
void
Propagator::resume(std::string const& checkpoint_directory, bool new_num_turns, int num_turns, bool new_max_turns,
        int max_turns, bool new_verbosity, int verbosity)
{
#if 0
    State state(get_resume_state(checkpoint_directory));
    if (new_max_turns) {
        state.max_turns = max_turns;
    }
    if (new_num_turns) {
    	state.num_turns = num_turns;
    }
    if (new_verbosity) {
        state.verbosity = verbosity;
    }
    if (state.bunch_simulator_ptr) {
        state.propagate_actions_ptr->before_resume_action(*stepper_sptr,
                                                          state.bunch_simulator_ptr->get_bunch());
    } else {
        state.propagate_actions_ptr->before_resume_action(*stepper_sptr,
                                                          state.bunch_train_simulator_ptr->get_bunch_train());
    }
    propagate(state);
    state.bunch_simulator_ptr ?  delete state.bunch_simulator_ptr: delete state.bunch_train_simulator_ptr;
    delete state.propagate_actions_ptr;
#endif
}

void
Propagator::propagate(
        Bunch_simulator & bunch_simulator, 
        int num_turns, int max_turns, int verbosity )
{
    Propagate_actions empty_propagate_actions;
    propagate( bunch_simulator, 
            empty_propagate_actions, 
            num_turns, max_turns, verbosity );
}

void
Propagator::propagate(
        Bunch_simulator & bunch_simulator,
        Propagate_actions & general_actions, 
        int num_turns, int max_turns, int verbosity )
{
    State state(&bunch_simulator, &general_actions, num_turns, 0, max_turns,
            verbosity);
    propagate(state);
}
#endif

#if 0
template<class Archive>
    void
    Propagator::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(stepper_sptr);
        ar & BOOST_SERIALIZATION_NVP(checkpoint_period);
        ar & BOOST_SERIALIZATION_NVP(checkpoint_dir);
        ar & BOOST_SERIALIZATION_NVP(checkpoint_with_xml);
        ar & BOOST_SERIALIZATION_NVP(concurrent_io);
        ar & BOOST_SERIALIZATION_NVP(final_checkpoint);
        ar & BOOST_SERIALIZATION_NVP(omp_threads);
    }

template
void
Propagator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Propagator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Propagator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Propagator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
#endif

