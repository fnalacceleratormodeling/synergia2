#include "standard_diagnostics_actions.h"

Standard_diagnostics_actions::Standard_diagnostics_actions()
{
}

void
Standard_diagnostics_actions::add_per_turn(Diagnostics_sptr diagnostics_sptr)
{
    per_turns.push_back(diagnostics_sptr);
}

void
Standard_diagnostics_actions::add_per_step(Diagnostics_sptr diagnostics_sptr)
{
    per_steps.push_back(diagnostics_sptr);
}

std::list<Diagnostics_sptr >
Standard_diagnostics_actions::get_per_steps_diagnostics_list() const
{
    return per_steps;
}

std::list<Diagnostics_sptr >
Standard_diagnostics_actions::get_per_turns_diagnostics_list() const
{
    return per_turns;
}

void
Standard_diagnostics_actions::update_and_write_all(
        std::list<Diagnostics_sptr > & diag_list)
{
    for (std::list<Diagnostics_sptr >::iterator it = diag_list.begin(); it
            != diag_list.end(); ++it) {
        (*it)->update_and_write();
    }
}

void
Standard_diagnostics_actions::first_action(Stepper & stepper, Bunch & bunch)
{
    update_and_write_all(per_turns);
    update_and_write_all(per_steps);
}

void
Standard_diagnostics_actions::turn_end_action(Stepper & stepper, Bunch & bunch,
        int turn_num)
{
    update_and_write_all(per_turns);
}

void
Standard_diagnostics_actions::step_end_action(Stepper & stepper, Step & step,
        Bunch & bunch, int turn_num, int step_num)
{
    update_and_write_all(per_steps);
}

Standard_diagnostics_actions::~Standard_diagnostics_actions()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Standard_diagnostics_actions)

