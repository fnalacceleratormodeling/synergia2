#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/independent_operation.h"


Bunch_simulator::Bunch_simulator(Bunch const & bunch)
    : num_turns(0)
    , first_turn(0)
    , max_turns(0)
    , trains()
    , diags_step()
    , diags_loss()
    , diags_ele()
    , diags_opr()
    , diags_opn()
{
}

void
Bunch_simulator::diag_action_step_and_turn(int turn_num, int step_num)
{
}

void
Bunch_simulator::diag_action_particle_loss_update()
{
}

void
Bunch_simulator::diag_action_particle_loss_write()
{
}

void
Bunch_simulator::diag_action_element(Lattice_element const & element)
{
}


void
Bunch_simulator::diag_action_operator(Operator const & opr)
{
}

void
Bunch_simulator::diag_action_operation(Independent_operation const & opn)
{
}

void
Bunch_simulator::prop_action_first(Lattice & lattice)
{
}

void
Bunch_simulator::prop_action_step_end(Lattice & lattice, int turn, int step)
{
}

void
Bunch_simulator::prop_action_turn_end(Lattice & lattice, int turn)
{
}

void 
Bunch_simulator::set_lattice_reference_particle(Reference_particle const & ref)
{
}
