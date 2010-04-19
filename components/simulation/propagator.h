#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "components/simulation/stepper.h"
#include "components/bunch/bunch.h"

class Propagator
{
private:
    Stepper_sptr stepper_sptr;

    void
    construct();
public:
    Propagator(Stepper_sptr const& stepper_sptr);
    void
    propagate(Bunch & bunch, int num_turns, bool diagnostics_per_step,
            bool diagnostics_per_turn);
    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
