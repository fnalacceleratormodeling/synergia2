#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "synergia/simulation/stepper.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics_writer.h"

class Propagator
{
private:
    Stepper_sptr stepper_sptr;

    void
    construct();
public:
    Propagator(Stepper_sptr stepper_sptr);
    void
    propagate(Bunch & bunch, int num_turns,
            Diagnostics & per_step_diagnostics,
            Diagnostics & per_turn_diagnostics, bool verbose = false);
    void
    propagate(Bunch & bunch, int num_turns,
            Multi_diagnostics & per_step_diagnostics,
            Multi_diagnostics & per_turn_diagnostics, bool verbose =
                    false);
    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
