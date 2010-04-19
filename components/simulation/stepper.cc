#include "stepper.h"
#include "utils/floating_point.h"

Steps &
Stepper::get_steps()
{
    return steps;
}

void
Stepper::print() const
{
    int index = 0;
    for (Steps::const_iterator it = steps.begin(); it != steps.end(); ++it) {
        ++index;
        (*it)->print(index);
    }
}

Stepper::~Stepper()
{

}

// Return an Independent_operator for a half step, starting at the
// lattice_element given by lattice_it at position left. Both lattice_it
// and left are updated by the function.
Independent_operator_sptr
Split_operator_stepper::get_half_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double half_step_length)
{
    Independent_operator_sptr retval(new Independent_operator(name,
            lattice_simulator.get_operation_extractor_map_sptr()));
    const double tolerance = 1.0e-8;
    double length = 0.0;
    bool complete = false;
    while (!complete) {
        double right = (*lattice_it)->get_length();
        if (floating_point_leq(length + (right - left), half_step_length,
                tolerance)) {
            Lattice_element_slice_sptr slice(new Lattice_element_slice(
                    *(*lattice_it), left, right));
            retval->append_slice(slice);
            length += (right - left);
            ++lattice_it;
            left = 0.0;
            if (floating_point_equal(length, half_step_length, tolerance)) {
                complete = true;
            } else {
                if (lattice_it == lattice_end) {
                    throw(std::runtime_error(
                            "get_half_step stepped beyond end of lattice"));
                }
            }
        } else {
            right = half_step_length - length + left;
            Lattice_element_slice_sptr slice(new Lattice_element_slice(
                    *(*lattice_it), left, right));
            retval->append_slice(slice);
            left = right;
            complete = true;
        }
    }
    return retval;
}

// extract_slices is a local function
Lattice_element_slices
extract_slices(Steps const& steps)
{
    Lattice_element_slices all_slices;
    for (Steps::const_iterator s_it = steps.begin(); s_it != steps.end(); ++s_it) {
        for (Operators::const_iterator o_it = (*s_it)->get_operators().begin(); o_it
                != (*s_it)->get_operators().end(); ++o_it) {
            if ((*o_it)->get_type() == "independent") {
                Lattice_element_slices
                        element_slices(boost::static_pointer_cast<
                                Independent_operator >(*o_it)->get_slices());
                all_slices.splice(all_slices.end(), element_slices);
            }
        }
    }
    return all_slices;
}

void
Split_operator_stepper::construct(
        Collective_operators const& collective_operators, int num_steps)
{
    double step_length = lattice_simulator.get_lattice_sptr()->get_length()
            / num_steps;
    double half_step_length = 0.5 * step_length;
    Lattice_elements::iterator lattice_it =
            lattice_simulator.get_lattice_sptr()->get_elements().begin();
    Lattice_elements::iterator lattice_end =
            lattice_simulator.get_lattice_sptr()->get_elements().end();
    double left = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        Step_sptr step(new Step);
        step->append(get_half_step("first_half", lattice_it, left, lattice_end,
                half_step_length));
        for (Collective_operators::const_iterator coll_op_it =
                collective_operators.begin(); coll_op_it
                != collective_operators.end(); ++coll_op_it) {
            step->append(*coll_op_it);
        }
        step->append(get_half_step("second_half", lattice_it, left,
                lattice_end, half_step_length));
        get_steps().push_back(step);
    }
    if (lattice_it != lattice_end) {
        throw(std::runtime_error(
                "internal error: split_operator_stepper did not make it to the end of the lattice\n"));
    }
    lattice_simulator.construct_sliced_chef_beamline(
            extract_slices(get_steps()));
}

Split_operator_stepper::Split_operator_stepper(
        Lattice_simulator const& lattice_simulator,
        Collective_operator_sptr collective_operator, int num_steps) :
    lattice_simulator(lattice_simulator)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, num_steps);
}

Split_operator_stepper::Split_operator_stepper(
        Lattice_simulator const& lattice_simulator,
        Collective_operators const& collective_operators, int num_steps) :
    lattice_simulator(lattice_simulator)
{
    construct(collective_operators, num_steps);
}

Split_operator_stepper::~Split_operator_stepper()
{

}


