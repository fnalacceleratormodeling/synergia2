#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagate_actions.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/train.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/logger.h"

class Propagator
{
public:
    static const char default_checkpoint_dir[];
    static const char description_file_name[];
    static const char propagator_archive_name[];
    static const char propagator_xml_archive_name[];
    static const char state_archive_name[];
    static const char state_xml_archive_name[];
    static const char log_file_name[];
    static const char stop_file_name[];
    static const char alt_stop_file_name[];

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
        {
        }
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version);
    };

private:
    Stepper_sptr stepper_sptr;
    int checkpoint_period;
    std::string checkpoint_dir;
    bool checkpoint_with_xml;

    void
    construct();
    State
    get_resume_state(std::string const& checkpoint_directory);
    void
    checkpoint(State & state, Logger & logger, double & t);
public:
    Propagator(Stepper_sptr stepper_sptr);

    // Default constructor for serialization use only
    Propagator();

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
    propagate(State & state);

    void
    resume(std::string const& checkpoint_dir, bool new_max_turns,
            int max_turns, bool new_verbosity, int verbosity);

    void
    propagate(Bunch_simulator & bunch_simulator, int num_turns,
            int max_turns = 0, int verbosity = 1);

    void
    propagate(Bunch_simulator & bunch_simulator,
            Propagate_actions & general_actions, int num_turns,
            int max_turns = 0, int verbosity = 1);

#if 0
    void
    propagate(Bunch_with_diagnostics_train & bunch_diag_train, int num_turns,
            bool verbose = false);

    void
    propagate(Bunch_with_diagnostics_train & bunch_diag_train, int num_turns,
            Propagate_actions & general_actions, bool verbose = true);
#endif

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
