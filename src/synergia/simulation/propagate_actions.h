#ifndef PROPAGATE_ACTIONS_H_
#define PROPAGATE_ACTIONS_H_

//#include "synergia/simulation/stepper.h"
#include "synergia/bunch/bunch.h"

class Stepper;
class Step;

class Propagate_actions
{
public:
    Propagate_actions();
    virtual void
    first_action(Stepper & stepper, Bunch & bunch);
    virtual void
    turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num);
    virtual void
    step_end_action(Stepper & stepper, Step & step, Bunch & bunch,
            int turn_num, int step_num);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
        }
    virtual
    ~Propagate_actions();
};

typedef boost::shared_ptr<Propagate_actions > Propagate_actions_sptr; // syndoc:include

#endif /* PROPAGATE_ACTIONS_H_ */
