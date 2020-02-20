#ifndef SYNERGIA_SIMULATION_CHECKPOINT_H
#define SYNERGIA_SIMULATION_CHECKPOINT_H

#include <utility>

class Propagator;
class Bunch_simulator;

namespace syn
{
    void checkpoint_save(Propagator const& prop, Bunch_simulator const& sim);

    std::pair<Propagator, Bunch_simulator>
    checkpoint_resume();
}

#endif
