#ifndef STANDARD_DIAGNOSTICS_ACTIONS_H_
#define STANDARD_DIAGNOSTICS_ACTIONS_H_

#include "synergia/simulation/propagate_actions.h"
#include "synergia/bunch/diagnostics.h"

class Standard_diagnostics_actions : public Propagate_actions
{
private:
    std::list<Diagnostics_sptr > per_turns, per_steps;
    void
    update_and_write_all(std::list<Diagnostics_sptr > & diag_list);
public:
    Standard_diagnostics_actions();
    void
    add_per_turn(Diagnostics_sptr diagnostics_sptr);
    void
    add_per_step(Diagnostics_sptr diagnostics_sptr);
    virtual void
    turn_start_action(Stepper & stepper, Bunch & bunch, int turn_num);
    virtual void
    step_start_action(Stepper & stepper, Step & step, Bunch & bunch,
            int turn_num, int step_num);
    virtual void
    final_action(Stepper & stepper, Bunch & bunch);
    virtual
    ~Standard_diagnostics_actions();

};
#endif /* STANDARD_DIAGNOSTICS_ACTIONS_H_ */
