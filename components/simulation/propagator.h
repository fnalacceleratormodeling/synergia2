#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "components/simulation/stepper.h"
#include "components/bunch/bunch.h"

class Propagator
{
private:
    Stepper stepper;

    void
    construct();
public:
    Propagator(Stepper & stepper);
    void
    propagate(Bunch & bunch, int num_turns, bool diagnostics_per_step,
            bool diagnostics_per_turn);
    ~Propagator();
};

#endif /* PROPAGATOR_H_ */
