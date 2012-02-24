#ifndef BUNCH_SIMULATOR_H_
#define BUNCH_SIMULATOR_H_

#include "synergia/bunch/bunch.h"
#include "synergia/simulation/standard_diagnostics_actions.h"

class Bunch_simulator
{
private:
    Bunch_sptr bunch_sptr;
    Standard_diagnostics_actions_sptr diagnostics_actions_sptr;

public:
    Bunch_simulator(Bunch_sptr bunch_sptr);
    Bunch_simulator(Bunch_sptr bunch_sptr,
            Standard_diagnostics_actions_sptr diagnostics_actions_sptr);
    // Default constructor for serialization use only
    Bunch_simulator();
    Bunch &
    get_bunch();
    Bunch_sptr
    get_bunch_sptr();
    Standard_diagnostics_actions &
    get_diagnostics_actions();
    Standard_diagnostics_actions_sptr
    get_diagnostics_actions_sptr();
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
