#include <iostream>

#include "step.h"

Step::Step() :
    operators()
{

}

void
Step::append(Operator_sptr operator_sptr)
{
    operators.push_back(operator_sptr);
}

void
Step::append(Operators const& the_operators)
{
    for (Operators::const_iterator it = the_operators.begin(); it
            != the_operators.end(); ++it) {
        operators.push_back(*it);
    }
}

Operators const&
Step::get_operators() const
{
    return operators;
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
