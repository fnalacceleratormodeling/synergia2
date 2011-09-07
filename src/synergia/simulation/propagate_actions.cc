#include "propagate_actions.h"

Propagate_actions::Propagate_actions()
{
}

void
Propagate_actions::turn_start_action(Stepper & stepper, Bunch & bunch,
        int turn_num)
{
}

void
Propagate_actions::step_start_action(Stepper & stepper, Step & step,
        Bunch & bunch, int turn_num, int step_num)
{
}

void
Propagate_actions::final_action(Stepper & stepper, Bunch & bunch)
{
}

Propagate_actions::~Propagate_actions()
{
}

Propagate_actions
empty_propagate_actions()
{
    return Propagate_actions();
}
