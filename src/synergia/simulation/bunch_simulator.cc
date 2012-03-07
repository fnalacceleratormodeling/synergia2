#include "bunch_simulator.h"
#include "diagnostics_actions.h"

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr) :
    bunch_sptr(bunch_sptr),
            diagnostics_actions_sptr(new Diagnostics_actions)
{
}

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr,
        Diagnostics_actions_sptr diagnostics_actions_sptr) :
    bunch_sptr(bunch_sptr), diagnostics_actions_sptr(diagnostics_actions_sptr)
{
}

Bunch_simulator::Bunch_simulator()
{
}

Bunch &
Bunch_simulator::get_bunch()
{
    return *bunch_sptr;
}

Bunch_sptr
Bunch_simulator::get_bunch_sptr()
{
    return bunch_sptr;
}

Diagnostics_actions &
Bunch_simulator::get_diagnostics_actions()
{
    return *diagnostics_actions_sptr;
}

Diagnostics_actions_sptr
Bunch_simulator::get_diagnostics_actions_sptr()
{
    return diagnostics_actions_sptr;
}

void
Bunch_simulator::add_per_turn(Diagnostics_sptr diagnostics_sptr, int period)
{
    diagnostics_actions_sptr->add_per_turn(diagnostics_sptr, period);
}

void
Bunch_simulator::add_per_turn(Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& turn_numbers)
{
    diagnostics_actions_sptr->add_per_turn(diagnostics_sptr, turn_numbers);
}

void
Bunch_simulator::add_per_step(Diagnostics_sptr diagnostics_sptr, int period)
{
    diagnostics_actions_sptr->add_per_step(diagnostics_sptr, period);
}

void
Bunch_simulator::add_per_step(Diagnostics_sptr diagnostics_sptr,
        std::list<int > const& step_numbers, int turn_period)
{
    diagnostics_actions_sptr->add_per_step(diagnostics_sptr, step_numbers,
            turn_period);
}

Bunch_simulator::~Bunch_simulator()
{
}
