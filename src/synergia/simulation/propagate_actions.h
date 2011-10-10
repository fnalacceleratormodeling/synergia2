#ifndef PROPAGATE_ACTIONS_H_
#define PROPAGATE_ACTIONS_H_

#include "synergia/simulation/stepper.h"
#include "synergia/bunch/bunch.h"

class Propagate_actions
{
public:
    Propagate_actions();
    virtual void
    first_action(Stepper & stepper, Bunch & bunch);
    virtual void
    turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num);
    virtual void
    step_end_action(Stepper & stepper, Step & step, Bunch & bunch,
            int turn_num, int step_num);
    virtual
    ~Propagate_actions();
};

#endif /* PROPAGATE_ACTIONS_H_ */
