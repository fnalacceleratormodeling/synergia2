#include "bunch_simulator.h"
#include "standard_diagnostics_actions.h"

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr) :
    bunch_sptr(bunch_sptr),
            diagnostics_actions_sptr(new Standard_diagnostics_actions)
{
}

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr,
        Standard_diagnostics_actions_sptr diagnostics_actions_sptr) :
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

Standard_diagnostics_actions &
Bunch_simulator::get_diagnostics_actions()
{
    return *diagnostics_actions_sptr;
}

Standard_diagnostics_actions_sptr
Bunch_simulator::get_diagnostics_actions_sptr()
{
    return diagnostics_actions_sptr;
}

Bunch_simulator::~Bunch_simulator()
{
}
