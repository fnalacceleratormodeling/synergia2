#ifndef SPLIT_OPERATOR_STEPPER_ELEMENTS_H_
#define SPLIT_OPERATOR_STEPPER_ELEMENTS_H_

#include "synergia/simulation/stepper.h"

/// The Split_operator_stepper_elements class generates a constant number of
/// split-operator steps per thick element. Thin elements are assigned
/// a single step each. One or more collective effects are included per
/// step.
class Split_operator_stepper_elements : public Stepper
{
private:

    //std::unique_ptr<Collective_operator> coll_opr;
    int steps_per_element;

public:

    explicit Split_operator_stepper_elements(
               int steps_per_element = 1 )
    : steps_per_element(steps_per_element)
    //, coll_opr(std::make_unique<COL_OP>(opr))
    { }

    std::vector<Step> apply(Lattice const & lattice) const override;
};

#endif /* SPLIT_OPERATOR_STEPPER_ELEMENTS_H_ */
