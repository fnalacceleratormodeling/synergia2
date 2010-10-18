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
            Diagnostics_writer & per_step_diagnostics,
            Diagnostics_writer & per_turn_diagnostics);
    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
