#ifndef SPLIT_OPERATOR_STEPPER_H_
#define SPLIT_OPERATOR_STEPPER_H_

#include "synergia/simulation/stepper.h"

/// The Split_operator_stepper class generates evenly-spaced split-operator
/// steps through a Lattice. One or more collective effects are included per
/// step.
class Split_operator_stepper : public Stepper
{
private:

    int num_steps;
    std::unique_ptr<CO_options> co_ops;

    std::vector<Step> 
        apply_impl(Lattice const & lattice) const override;


public:

    Split_operator_stepper(CO_options const & coo, int num_steps)
    : num_steps(num_steps), co_ops(coo.clone())
    { }
};


#endif /* SPLIT_OPERATOR_STEPPER_H_ */
