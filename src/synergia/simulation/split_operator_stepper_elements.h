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

    int steps_per_element;
    std::unique_ptr<CO_options> co_ops;

    std::vector<Step> 
        apply_impl(Lattice const & lattice) const override;

public:

    template<class COO>
    Split_operator_stepper_elements(
               int steps_per_element,
               COO const & coo )
    : steps_per_element(steps_per_element)
    , co_ops(std::make_unique<COO>(coo))
    { }
};

#endif /* SPLIT_OPERATOR_STEPPER_ELEMENTS_H_ */
