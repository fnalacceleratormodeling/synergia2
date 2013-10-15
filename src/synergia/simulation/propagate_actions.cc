#include "propagate_actions.h"

Propagate_actions::Propagate_actions()
{
}
Propagate_actions::Propagate_actions(std::string const& type):
type(type)
{
}  
void
Propagate_actions::first_action(Stepper & stepper, Bunch & bunch)
{
}

void
Propagate_actions::first_action(Stepper & stepper, Bunch_train & bunch_train)
{
}

void
Propagate_actions::turn_end_action(Stepper & stepper, Bunch & bunch,
        int turn_num)
{
}

void
Propagate_actions::turn_end_action(Stepper & stepper, Bunch_train & bunch_train,
        int turn_num)
{
}

void
Propagate_actions::step_end_action(Stepper & stepper, Step & step,
        Bunch & bunch, int turn_num, int step_num)
{
}

void
Propagate_actions::step_end_action(Stepper & stepper, Step & step,
        Bunch_train & bunch_train, int turn_num, int step_num)
{
}

void
Propagate_actions::operator_begin_action(Stepper & stepper, Step & step, Operator & op, int step_num, int turn_num, 
                   int bunch_num)                      
{
}  

std::string const&
Propagate_actions::get_type() const
{
  return type;
}  



Propagate_actions::~Propagate_actions()
{
}
