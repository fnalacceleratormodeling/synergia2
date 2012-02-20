#ifndef DIAGNOSTICS_ACTIONS_H_
#define DIAGNOSTICS_ACTIONS_H_

#include "synergia/simulation/stepper.h"
#include "synergia/bunch/bunch.h"

//class Stepper;
//class Step;

class Diagnostic_actions
{
public:
    Diagnostic_actions()
    {
    }
    virtual void
    first_action(Stepper & stepper, Bunch & bunch) = 0;
    virtual void
    turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num) = 0;
    virtual void
    step_end_action(Stepper & stepper, Step & step, Bunch & bunch,
            int turn_num, int step_num) = 0;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
        }
    virtual
    ~Diagnostic_actions()
    {
    }
};

typedef boost::shared_ptr<Diagnostic_actions > Diagnostic_actions_sptr;

#endif /* DIAGNOSTICS_ACTIONS_H_ */
