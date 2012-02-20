#include "bunch_simulator.h"
#include "standard_diagnostics_actions.h"

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr) :
    bunch_sptr(bunch_sptr),
            diagnostic_actions_sptr(new Standard_diagnostics_actions)
{
}

Bunch_simulator::Bunch_simulator(Bunch_sptr bunch_sptr,
        Diagnostic_actions_sptr diagnostic_actions_sptr) :
    bunch_sptr(bunch_sptr), diagnostic_actions_sptr(diagnostic_actions_sptr)
{
}

Bunch_sptr
Bunch_simulator::get_bunch_sptr()
{
    return bunch_sptr;
}

Diagnostic_actions_sptr
Bunch_simulator::get_diagnostic_actions_sptr()
{
    return diagnostic_actions_sptr;
}

Bunch_simulator::~Bunch_simulator()
{
}
