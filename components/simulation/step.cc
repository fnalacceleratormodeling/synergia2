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
    Operators tmp(the_operators);
    operators.splice(operators.end(), tmp);
}

void
Step::apply(Bunch & bunch)
{
    for (Operators::const_iterator it = operators.begin(); it
            != operators.end(); ++it) {
        (*it)->apply(bunch, operators);
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
