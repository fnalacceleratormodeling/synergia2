#include <iostream>

#include "step.h"
#include "synergia/foundation/physical_constants.h"

Step::Step(double length) :
    operators(), time_fractions(), length(length)
{

}

void
Step::append(Operator_sptr operator_sptr, double time_fraction)
{
    operators.push_back(operator_sptr);
    time_fractions.push_back(time_fraction);
}

void
Step::append(Operators const& the_operators, double time_fraction)
{
    Operators tmp(the_operators);
    operators.splice(operators.end(), tmp);
    for (Operators::const_iterator it = the_operators.begin(); it
            != the_operators.end(); ++it) {
        time_fractions.push_back(time_fraction);
    }
}

void
Step::apply(Bunch & bunch)
{
    std::list<double >::const_iterator fractions_it = time_fractions.begin();
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        // time [s] in accelerator frame
        double time = length / (bunch.get_reference_particle().get_beta()
                * pconstants::c);
        (*it)->apply(bunch, (*fractions_it) * time, *this);
        ++fractions_it;
    }
}

Operators const&
Step::get_operators() const
{
    return operators;
}

std::list<double > const&
Step::get_time_fractions() const
{
    return time_fractions;
}

void
Step::print(int index) const
{
    std::cout << "step " << index << ":\n";
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        (*it)->print();
    }
}
