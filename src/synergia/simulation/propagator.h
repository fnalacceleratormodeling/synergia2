#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagate_actions.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_with_diagnostics.h"
#include "synergia/bunch/train.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/utils/serialization.h"

class Propagator
{
private:
    Stepper_sptr stepper_sptr;

    void
    construct();
public:
    struct State
    {
        Bunch_with_diagnostics * bunch_with_diagnostics_ptr;
        int num_turns;
        int first_turn;
        Propagate_actions * general_actions_ptr;
        bool verbose;
        State(Bunch_with_diagnostics * bunch_with_diagnostics_ptr, int num_turns,
                int first_turn, Propagate_actions * propagate_actions_ptr, bool verbose);
        State(){}
        template<class Archive>
            void
            serialize(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(bunch_with_diagnostics_ptr);
                ar & BOOST_SERIALIZATION_NVP(num_turns);
                ar & BOOST_SERIALIZATION_NVP(first_turn);
                ar & BOOST_SERIALIZATION_NVP(general_actions_ptr);
                ar & BOOST_SERIALIZATION_NVP(verbose);
            }
    };
    Propagator(Stepper_sptr stepper_sptr);

    // Default propagator for serialization use only
    Propagator();

    void
    propagate(State & state);

    void
    propagate(Bunch_with_diagnostics & bunch_with_diagnostics, int num_turns,
            bool verbose = false);

    void
    propagate(Bunch_with_diagnostics & bunch_with_diagnostics, int num_turns,
            Propagate_actions & general_actions, bool verbose = false);

    void
    propagate(Bunch_with_diagnostics_train & bunch_diag_train, int num_turns,
            bool verbose = false);

    void
    propagate(Bunch_with_diagnostics_train & bunch_diag_train, int num_turns,
            Propagate_actions & general_actions, bool verbose = false);

    void
    propagate(Bunch & bunch, int num_turns,
            Diagnostics & per_step_diagnostics,
            Diagnostics & per_turn_diagnostics,
            bool verbose = false);

    void
    propagate(Bunch & bunch, int num_turns,
            Multi_diagnostics & per_step_diagnostics,
            Multi_diagnostics & per_turn_diagnostics, bool verbose = false);

    void
    propagate(Bunch & bunch, int num_turns,
            Standard_diagnostics_actions & diagnostics_actions,
            int verbosity = 0);

    void
    propagate(Bunch & bunch, int num_turns,
            Standard_diagnostics_actions & diagnostics_actions,
            Propagate_actions & general_actions, int verbosity = 0);

     template<class Archive>
         void
         serialize(Archive & ar, const unsigned int version)
         {
             ar & BOOST_SERIALIZATION_NVP(stepper_sptr);
         }


    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
