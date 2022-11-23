#include "synergia/simulation/propagator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/checkpoint.h"

#include "synergia/utils/digits.h"

void
Propagator::do_before_start(Bunch_simulator& simulator, Logger& logger)
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

    if (simulator.current_turn() == 0) {
        Kokkos::Profiling::pushRegion("diagnostic-actions");
        simulator.diag_action_step_and_turn(PRE_TURN, FINAL_STEP);
        Kokkos::Profiling::popRegion();

        simulator.prop_action_first(lattice);
        simulator.set_lattice_reference_particle(
            lattice.get_reference_particle());
    }

    auto updates = lattice.update();
    // if (updates.structure) steps = stepper.apply(lattice);
}

void
Propagator::do_start_repetition(Bunch_simulator& simulator)
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

    for (auto& bunch : simulator[0].get_bunches())
        bunch.get_reference_particle().start_repetition();

    for (auto& bunch : simulator[1].get_bunches())
        bunch.get_reference_particle().start_repetition();
}

void
Propagator::do_step(Bunch_simulator& simulator,
                    Step& step,
                    int step_count,
                    int turn_count,
                    Logger& logger)
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
    Kokkos::Profiling::pushRegion("diagnostic-actions");
    simulator.diag_action_step_and_turn(turn_count, step_count);
    Kokkos::Profiling::popRegion();

    // t = simple_timer_show(t, "propagate-diagnostics_actions_step");

    // propagate action
    simulator.prop_action_step_end(lattice, turn_count, step_count);
    // t = simple_timer_show(t, "propagate-general_actions-step");

    double t_step1 = MPI_Wtime();

    logger(LoggerV::INFO_STEP)
        << "Propagator: step " << std::setw(digits(steps.size())) << step_count
        << "/" << steps.size()

        << ", s_n = " << std::fixed << std::setprecision(4)
        << simulator[0][0].get_reference_particle().get_s_n()

        << ", time = " << std::fixed << std::setprecision(3)
        << t_step1 - t_step0 << "s, macroparticles = ";

    for (auto const& train : simulator.get_trains()) {
        logger << "(";

        for (auto const& bunch : train.get_bunches()) {
            logger << bunch.get_total_num();
            if (bunch.get_array_index() != train.get_bunch_array_size() - 1)
                logger << ", ";
        }

        logger << ")";
        if (train.get_index() == 0) logger << " / ";
    }

    logger << "\n";

    logger(LoggerV::INFO_OPR) << "\n\n";

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
Propagator::check_out_of_particles(Bunch_simulator const& simulator,
                                   Logger& logger)
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
Propagator::do_turn_end(Bunch_simulator& simulator,
                        int turn_count,
                        Logger& logger)
{
    // t = simple_timer_current();

    // diagnostic actions
    Kokkos::Profiling::pushRegion("diagnostic-actions");
    simulator.diag_action_step_and_turn(turn_count, Propagator::FINAL_STEP);
    Kokkos::Profiling::popRegion();

    // propagate actions
    simulator.prop_action_turn_end(lattice, turn_count);

    // update lattice in case it has been chaged in the propagate action
    auto updates = lattice.update();
    // if (updates.structure()) steps = stepper.apply(lattice);

    // increment the turn number
    simulator.inc_turn();

    // double t_turn1 = MPI_Wtime();

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
Propagator::propagate(Bunch_simulator& sim, Logger& logger, int max_turns)
{
    Kokkos::Profiling::pushRegion("propagate-execute");
    const int total_turns = sim.max_turns();

    // parameter check
    if (max_turns == -1 && total_turns == -1) {
        throw std::runtime_error(
            "Max number of turns must be set either in the Bunch_simulator "
            "(Bunch_simulator::set_num_turns(int)), or in the Propagator "
            "(Propagator::propagate(..., int max_turns))");
    }

    try {
        do_before_start(sim, logger);

        // first turn is always the current turn from simulator
        int turn = sim.current_turn();

        // decide the last turn
        int last_turn = (max_turns == -1) ? total_turns : turn + max_turns;
        if ((last_turn > total_turns) && (total_turns != -1))
            last_turn = total_turns;

        int turns_since_checkpoint = 0;
        bool out_of_particles = false;
        double t_prop0 = MPI_Wtime();

        logger(LoggerV::INFO_TURN) << "Propagator: starting turn " << turn + 1
                                   << ", final turn " << last_turn << "\n\n";

        for (; turn < last_turn; ++turn) {
            double t_turn0 = MPI_Wtime();

            do_start_repetition(sim);

            int step_count = 0;
            for (auto& step : steps) {
                ++step_count;
                do_step(sim, step, step_count, turn, logger);

                out_of_particles = check_out_of_particles(sim, logger);
                if (out_of_particles) break;
            }

            double t_turn1 = MPI_Wtime();

            // turn end log
            logger(LoggerV::INFO_TURN)
                << "Propagator: turn " << std::setw(digits(total_turns))
                << turn + 1 << "/";

            if (total_turns == -1)
                logger << "inf.";
            else
                logger << total_turns;

            logger << ", time = " << std::fixed << std::setprecision(3)
                   << t_turn1 - t_turn0 << "s"

                   << ", macroparticles = ";

            for (auto const& train : sim.get_trains()) {
                logger << "(";

                for (auto const& bunch : train.get_bunches()) {
                    logger << bunch.get_total_num();
                    if (bunch.get_array_index() !=
                        train.get_bunch_array_size() - 1)
                        logger << ", ";
                }

                logger << ")";
                if (train.get_index() == 0) logger << " / ";
            }

            logger << "\n";
            logger(LoggerV::INFO_STEP) << "\n";

            // out of particles
            if (out_of_particles) break;

            ++turns_since_checkpoint;
            do_turn_end(sim, turn, logger);

            // checkpoint save
            // syn::checkpoint_save(*this, sim);

            if ((turns_since_checkpoint == checkpoint_period) ||
                ((turn == (sim.max_turns() - 1)) && final_checkpoint)) {
                // t = simple_timer_current();
                syn::checkpoint_save(*this, sim);
                // t = simple_timer_show(t, "propagate-checkpoint_period");
                turns_since_checkpoint = 0;
            }

#if 0
            if (((turn - orig_first_turn + 1) == sim.max_turns) 
                    && (turn != (sim.num_turns - 1))) 
            {
                logger(LoggerV::INFO_TURN) 
                    << "Propagator: maximum number of turns reached\n";

                if (turns_since_checkpoint > 0) 
                {
                    t = simple_timer_current();
                    checkpoint(sim, logger, t);
                    t = simple_timer_show(t, "propagate-checkpoint_max");
                }

                break;
            }
#endif

#if 0
            if (fs::exists(stop_file_name) || fs::exists(alt_stop_file_name)) 
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

        if (last_turn != total_turns) {
            logger(LoggerV::INFO_TURN)
                << "Propagator: maximum number of turns reached\n";

            // TODO: checkpoint
        }

        if (out_of_particles) {
            logger(LoggerV::WARNING) << "Propagator: no particles left\n";
        }

        double t_prop1 = MPI_Wtime();

        logger(LoggerV::INFO_TURN)
            << "Propagator: total time = " << std::fixed << std::setprecision(3)
            << t_prop1 - t_prop0 << "s\n";
    }
    catch (std::exception const& e) {
        std::cerr << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 888);
    }

    Kokkos::Profiling::popRegion();
}
