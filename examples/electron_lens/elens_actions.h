#ifndef __ELENS_ACTIONS__
#define __ELENS_ACTIONS__

#include "synergia/simulation/propagate_actions.h"

class Elens_actions: public Propagate_actions {
public:
    Elens_actions(); // for serialization use
    Elens_actions(bool elensadaptive);
    void first_action(Stepper & stepper, Bunch & bunch);
    void turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num);
    void first_action(Stepper & stepper, Bunch_train & bunch_train);
    void turn_end_action(Stepper & stepper, Bunch_train & bunch_train, int turn_num);
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version);
    ~Elens_actions();
private:
    bool adaptive;
    void set_elens_adaptive_radius(Stepper & stepper, Bunch & bunch);
    void set_elens_adaptive_radius(Stepper & stepper, Bunch_train & bunch_train);
};

BOOST_CLASS_EXPORT_KEY(Elens_actions)

#endif // __ELENS_ACTIONS__
