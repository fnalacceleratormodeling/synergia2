#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagate_actions.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/train.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/utils/serialization.h"

class Propagator
{
public:
    static const char default_checkpoint_dir[];
    static const char propagator_archive_name[];
    static const char state_archive_name[];

    struct State
    {
        Bunch_simulator * bunch_simulator_ptr;
        Propagate_actions * propagate_actions_ptr;
        int num_turns;
        int first_turn;
        int max_turns;
        bool verbose;
        State(Bunch_simulator * bunch_simulator_ptr,
                Propagate_actions * propagate_actions_ptr, int num_turns,
                int first_turn, int max_turns, bool verbose);
        State()
        {
        }
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(bunch_simulator_ptr);
                ar & BOOST_SERIALIZATION_NVP(propagate_actions_ptr);
                ar & BOOST_SERIALIZATION_NVP(num_turns);
                ar & BOOST_SERIALIZATION_NVP(first_turn);
                ar & BOOST_SERIALIZATION_NVP(max_turns);
                ar & BOOST_SERIALIZATION_NVP(verbose);
            }
    };

private:
    Stepper_sptr stepper_sptr;
    int checkpoint_period;
    std::string checkpoint_dir;

    void
    construct();
    State
    get_resume_state(std::string const& checkpoint_directory);

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
    propagate(State & state);

    void
    checkpoint(State & state);

    void
    resume(std::string const& checkpoint_dir);

    void
    resume(std::string const& checkpoint_dir, int max_turns);

    void
    propagate(Bunch_simulator & bunch_simulator, int num_turns,
            int max_turns = 0, bool verbose = true);

    void
    propagate(Bunch_simulator & bunch_simulator,
            Propagate_actions & general_actions, int num_turns,
            int max_turns = 0, bool verbose = true);

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
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(stepper_sptr);
            ar & BOOST_SERIALIZATION_NVP(checkpoint_period);
            ar & BOOST_SERIALIZATION_NVP(checkpoint_dir);
        }

    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
