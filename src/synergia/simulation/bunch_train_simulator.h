#ifndef BUNCH_TRAIN_SIMULATOR_H_
#define BUNCH_TRAIN_SIMULATOR_H_

#include "synergia/bunch/bunch_train.h"
#include "synergia/simulation/diagnostics_actions.h"

class Bunch_train_simulator
{
private:
    Bunch_train_sptr bunch_train_sptr;
    Diagnostics_actionss diagnostics_actionss;

public:
    Bunch_train_simulator(Bunch_train_sptr bunch_train_sptr);
    // Default constructor for serialization use only
    Bunch_train_simulator();
    Bunch_train &
    get_bunch_train();
    Bunch_train_sptr
    get_bunch_train_sptr();
    Diagnostics_actionss &
    get_diagnostics_actionss();
    void
    add_per_turn(int which, Diagnostics_sptr diagnostics_sptr, int period = 1);
    void
    add_per_turn(int which, Diagnostics_sptr diagnostics_sptr,
            std::list<int > const& turn_numbers);
    void
    add_per_step(int which, Diagnostics_sptr diagnostics_sptr, int period = 1);
    void
    add_per_step(int which, Diagnostics_sptr diagnostics_sptr,
            std::list<int > const& step_numbers, int turn_period = 1);
    void
    add_per_forced_diagnostics_step(int which,
            Diagnostics_sptr diagnostics_sptr, int turn_period = 1);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Bunch_train_simulator();
};

#endif /* BUNCH_TRAIN_SIMULATOR_H_ */
