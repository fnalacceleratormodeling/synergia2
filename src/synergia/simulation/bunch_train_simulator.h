#ifndef BUNCH_TRAIN_SIMULATOR_H_
#define BUNCH_TRAIN_SIMULATOR_H_

#include "synergia/bunch/bunch.h"
#include "synergia/simulation/diagnostics_actions.h"

class Bunch_train_simulator
{
private:
    Bunches bunches;
    std::vector<double > spacings;
    Diagnostics_actionss diagnostics_actionss;

public:
    Bunch_train_simulator(Bunches const& bunches, double spacing);
    Bunch_train_simulator(Bunches const& bunches,
            std::vector<double > const& spacings);
    // Default constructor for serialization use only
    Bunch_train_simulator();
    size_t
    get_size() const;
    Bunches &
    get_bunches();
    std::vector<double > &
    get_spacings();
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
