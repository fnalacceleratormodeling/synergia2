#ifndef INDEPENDENT_STEPPER_ELEMENTS_H_
#define INDEPENDENT_STEPPER_ELEMENTS_H_

#include "synergia/simulation/stepper.h"

/// The Independent_stepper_elements class generates a constant number of
/// Independent_operator steps per thick element. Thin elements are assigned
/// a single step each. No collective effects are included.
class Independent_stepper_elements : public Stepper
{
private:

    int steps_per_element;

    std::vector<Step> 
    apply_impl(Lattice const & lattice) const override;

public:

    Independent_stepper_elements(int steps_per_element)
        : steps_per_element(steps_per_element)
    { }
};

#endif /* INDEPENDENT_STEPPER_ELEMENTS_H_ */
