#include "stepper.h"
#include "utils/floating_point.h"

Independent_operator::Independent_operator() :
    slice_list()
{

}

void
Independent_operator::append_slice(Lattice_element_slice const& slice)
{
    slice_list.push_back(slice);
}

Lattice_element_slices const&
Independent_operator::get_slices() const
{
    return slice_list;
}

void
Independent_operator::print() const
{
    std::cout << "Independent_operator:\n";
}

Independent_operator::~Independent_operator()
{

}

Step::Step() :
    operators()
{

}

void
Step::append(Operator const& the_operator)
{
    operators.push_back(the_operator);
}

void
Step::append(Operators const& the_operators)
{
    for (Operators::const_iterator it = the_operators.begin(); it
            != the_operators.end(); ++it) {
        operators.push_back(*it);
    }
}

void
Step::print() const
{
    std::cout << "step:\n";
    for(Operators::const_iterator it = operators.begin(); it != operators.end(); ++it) {
        it->print();
    }
}

Steps &
Stepper::get_steps()
{
    return steps;
}

Lattice &
Stepper::get_lattice()
{
    return *lattice_ptr;
}

Chef_lattice &
Stepper::get_chef_lattice()
{
    return *chef_lattice_ptr;
}

// Return an Independent_operator for a half step, starting at the
// lattice_element given by lattice_it at position left. Both lattice_it
// and left are updated by the function.
Independent_operator
Split_operator_stepper::get_half_step(
        Lattice_element_list::iterator & lattice_it, double & left,
        Lattice_element_list::iterator const & lattice_end,
        const double half_step_length)
{
    Independent_operator retval;
    const double tolerance = 1.0e-8;
    double length = 0.0;
    bool complete = false;
    while (!complete) {
        double right = lattice_it->get_length();
        if (floating_point_leq(length + (right - left), half_step_length,
                tolerance)) {
            Lattice_element_slice slice(*lattice_it, left, right);
            retval.append_slice(slice);
            length += (right - left);
            ++lattice_it;
            if (floating_point_equal(length, half_step_length, tolerance)) {
                complete = true;
            } else {
                if (lattice_it == lattice_end) {
                    throw(std::runtime_error(
                            "get_half_step stepped beyond end of lattice"));
                }
                left = 0.0;
            }
        } else {
            right = half_step_length - length + left;
            Lattice_element_slice slice(*lattice_it, left, right);
            retval.append_slice(slice);
            left = right;
            complete = true;
        }
    }
    return retval;
}

void
Split_operator_stepper::construct(Lattice & lattice,
        Chef_lattice & chef_lattice, int num_steps,
        Collective_operators const& collective_operators)
{
    lattice_ptr = &lattice;
    chef_lattice_ptr = &chef_lattice;

    double step_length = lattice_ptr->get_length() / num_steps;
    double half_step_length = 0.5 * step_length;
    Lattice_element_list::iterator lattice_it = lattice.get_elements().begin();
    Lattice_element_list::iterator lattice_end = lattice.get_elements().end();
    double left = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        Step step;
        step.append(get_half_step(lattice_it, left, lattice_end,
                half_step_length));
        for (Collective_operators::const_iterator coll_op_it =
                collective_operators.begin(); coll_op_it
                != collective_operators.end(); ++coll_op_it) {
            step.append(*coll_op_it);
        }
        step.append(get_half_step(lattice_it, left, lattice_end,
                half_step_length));
        steps.push_back(step);
    }
    if (lattice_it == lattice_end) {
        std::cout
                << "jfa: all is right with the world, lattice_it == lattice_end after all split operator steps have been constructed\n";
    } else {
        std::cout
                << "jfa: sigh... split_operator_stepper didn't make it to the end of the lattice\n";
    }
}

Split_operator_stepper::Split_operator_stepper(Lattice & lattice,
        Chef_lattice & chef_lattice, int num_steps,
        Collective_operator const& collective_operator)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(lattice, chef_lattice, num_steps, collective_operators);
}

Split_operator_stepper::Split_operator_stepper(Lattice & lattice,
        Chef_lattice & chef_lattice, int num_steps,
        Collective_operators const& collective_operators)
{
    construct(lattice, chef_lattice, num_steps, collective_operators);
}

