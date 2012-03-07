#ifndef BUNCH_SIMULATOR_H_
#define BUNCH_SIMULATOR_H_

#include "synergia/bunch/bunch.h"
#include "synergia/simulation/diagnostics_actions.h"

class Bunch_simulator
{
private:
    Bunch_sptr bunch_sptr;
    Diagnostics_actions_sptr diagnostics_actions_sptr;

public:
    Bunch_simulator(Bunch_sptr bunch_sptr);
    Bunch_simulator(Bunch_sptr bunch_sptr,
            Diagnostics_actions_sptr diagnostics_actions_sptr);
    // Default constructor for serialization use only
    Bunch_simulator();
    Bunch &
    get_bunch();
    Bunch_sptr
    get_bunch_sptr();
    Diagnostics_actions &
    get_diagnostics_actions();
    Diagnostics_actions_sptr
    get_diagnostics_actions_sptr();
    void
    add_per_turn(Diagnostics_sptr diagnostics_sptr, int period = 1);
    void
    add_per_turn(Diagnostics_sptr diagnostics_sptr,
            std::list<int > const& turn_numbers);
    void
    add_per_step(Diagnostics_sptr diagnostics_sptr, int period = 1);
    void
    add_per_step(Diagnostics_sptr diagnostics_sptr,
            std::list<int > const& step_numbers, int turn_period = 1);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
            ar & BOOST_SERIALIZATION_NVP(diagnostics_actions_sptr);
        }
    ~Bunch_simulator();
};
#endif /* BUNCH_SIMULATOR_H_ */
