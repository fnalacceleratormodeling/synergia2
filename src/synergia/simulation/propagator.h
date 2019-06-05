#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/lattice/lattice.h"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"

#include "synergia/utils/cereal.h"
#include "synergia/utils/logger.h"

class Propagator
{

#if 0
public:
    static const std::string default_checkpoint_dir;
    static const std::string description_file_name;
    static const std::string propagator_archive_name;
    static const std::string propagator_xml_archive_name;
    static const std::string state_archive_name;
    static const std::string state_xml_archive_name;
    static const std::string log_file_name;
    static const std::string stop_file_name;
    static const std::string alt_stop_file_name;
    static const int default_checkpoint_period;
    static const int default_concurrent_io;

    struct State
    {
        Bunch_simulator * bunch_simulator_ptr;
        Propagate_actions * propagate_actions_ptr;

        int num_turns;
        int first_turn;
        int max_turns;
        int verbosity;

        State(Bunch_simulator * bunch_simulator_ptr,
                Propagate_actions * propagate_actions_ptr, int num_turns,
                int first_turn, int max_turns, int verbosity);

        State() 
        { }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    };

private:

    Stepper_sptr stepper_sptr;

    int checkpoint_period;
    std::string checkpoint_dir;
    bool checkpoint_with_xml;
    int concurrent_io;
    bool final_checkpoint;

    int omp_threads;
#endif

public:

    static const int FINAL_STEP = -1;

private:

    Lattice lattice;
    std::vector<Step> steps;

private:

	void construct();

	void do_before_start(
            Bunch_simulator & simulator, 
            Logger & logger);

	void do_start_repetition(
            Bunch_simulator & simulator);

	void do_turn_end(
            Bunch_simulator & simulator, 
            int turn_count, 
            Logger & logger);

	void do_step(
            Bunch_simulator & simulator,
            Step & step, 
            int step_count, 
            int turn_count, 
            Logger & logger);

	bool check_out_of_particles(
            Bunch_simulator const & simulator, 
            Logger & logger);

#if 0
	void
	checkpoint(State & state, Logger & logger, double & t);
#endif

public:

    explicit Propagator(
            Lattice const & lattice, 
            Stepper const & stepper = Split_operator_stepper_elements(1) );

    void propagate(Bunch_simulator & simulator);

    void print_steps(Logger & logger) const
    { for(auto const & s : steps) s.print(logger); }

#if 0
    Stepper_sptr
    get_stepper_sptr();

    void
    set_checkpoint_period(int period);

    int
    get_checkpoint_period() const;

    void
    set_checkpoint_dir(std::string const& directory_name);

    std::string const&
    get_checkpoint_dir() const;

    void
    set_checkpoint_with_xml(bool with_xml);

    bool
    get_checkpoint_with_xml() const;

    void
    set_final_checkpoint(bool final_checkpoint);

    bool
    get_final_checkpoint() const;

    void
    set_concurrent_io(int max);

    int
    get_concurrent_io() const;

    void
    set_num_threads(int nt);

    int
    get_num_threads() const;

    void
    propagate(State & state);

    /// jfa note: the lifetime of the pointers in state must
    ///           be managed manually
    State
    get_resume_state(std::string const& checkpoint_dir);

    void
    resume( std::string const& checkpoint_dir, 
            bool new_num_turns, int num_turns, 
            bool new_max_turns, int max_turns, 
            bool new_verbosity, int verbosity );

    void
    propagate(
            Bunch_simulator & bunch_simulator, 
            int num_turns,
            int max_turns = 0, 
            int verbosity = 1 );

    void
    propagate(
            Bunch_simulator & bunch_simulator,
            Propagate_actions & general_actions, 
            int num_turns,
            int max_turns = 0, 
            int verbosity = 1 );
#endif

#if 0
    template<class Archive>
    void
    serialize(Archive & ar, const unsigned int version);
#endif
};

#endif /* PROPAGATOR_H_ */
