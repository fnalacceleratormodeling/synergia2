#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "components/simulation/stepper.h"
#include "components/bunch/bunch.h"
#include "components/bunch/diagnostics_writer.h"

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
            Diagnostics_writer & per_step_diagnostics,
            Diagnostics_writer & per_turn_diagnostics);
    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
