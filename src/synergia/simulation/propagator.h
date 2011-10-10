#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagate_actions.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_with_diagnostics.h"
#include "synergia/bunch/train.h"
#include "synergia/bunch/multi_diagnostics.h"

class Propagator
{
private:
    Stepper_sptr stepper_sptr;

    void
    construct();
public:
    Propagator(Stepper_sptr stepper_sptr);
    
    void
    propagate(Bunch_with_diagnostics & bunch_with_diagnostics, int num_turns, 
              bool verbose = false);
     void
     propagate(Bunch_with_diagnostics_train & bunch_diag_train, int num_turns,
              bool verbose= false);
    
    void
    propagate(Bunch & bunch, int num_turns, Diagnostics & per_step_diagnostics,
            Diagnostics & per_turn_diagnostics, bool verbose = false);
    void
    propagate(Bunch & bunch, int num_turns,
            Multi_diagnostics & per_step_diagnostics,
            Multi_diagnostics & per_turn_diagnostics, bool verbose =
                    false);
    
    void
    propagate(Bunch & bunch, int num_turns,
            Propagate_actions & diagnostics_actions, int verbosity = 0);
    void
    propagate(Bunch & bunch, int num_turns,
            Propagate_actions & diagnostics_actions,
            Propagate_actions & general_actions, int verbosity = 0);

    
                
    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
