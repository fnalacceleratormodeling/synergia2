#include "stepper.h"
#include "synergia/utils/floating_point.h"

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

Independent_operator_sptr
Independent_stepper::get_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double step_length)
{
    Independent_operator_sptr retval(new Independent_operator(name,
            lattice_simulator.get_operation_extractor_map_sptr()));
    const double tolerance = 1.0e-8;
    double length = 0.0;
    bool complete = false;
    while (!complete) {
        double right = (*lattice_it)->get_length();
        if (floating_point_leq(length + (right - left), step_length,
                tolerance)) {
            // The rest of the element fits in the half step
            Lattice_element_slice_sptr slice(new Lattice_element_slice(
                    *(*lattice_it), left, right));
            retval->append_slice(slice);
            length += (right - left);
            ++lattice_it;
            left = 0.0;
            if (floating_point_equal(length, step_length, tolerance)) {
                if ((lattice_it == lattice_end) || ((*lattice_it)->get_length()
                        != 0.0)) {
                    complete = true;
                }
            } else {
                if (lattice_it == lattice_end) {
                    throw(std::runtime_error(
                            "get_step stepped beyond end of lattice"));
                }
            }
        } else {
            // Need to take a portion of the element...
            bool end_within_error = false;
            double old_right = right;
            right = step_length - length + left;
            if ((old_right - right) < retval->get_slices().size() * tolerance) {
                // ... unless we are within an accumulated tolerance of the end
                right = old_right;
                end_within_error = true;
            }
            Lattice_element_slice_sptr slice(new Lattice_element_slice(
                    *(*lattice_it), left, right));
            retval->append_slice(slice);
            length += (right - left);
            if (end_within_error) {
                ++lattice_it;
                left = 0.0;
            } else {
                left = right;
            }
            complete = true;
        }
    }
    return retval;
}

Independent_stepper::Independent_stepper(
        Lattice_simulator const& lattice_simulator, int num_steps) :
    lattice_simulator(lattice_simulator)
{
    double step_length = this->lattice_simulator.get_lattice_sptr()->get_length()
            / num_steps;

    Lattice_elements::iterator lattice_it =
            this->lattice_simulator.get_lattice_sptr()->get_elements().begin();
    Lattice_elements::iterator lattice_end =
            this->lattice_simulator.get_lattice_sptr()->get_elements().end();

    double left = 0.0;

    for (int i = 0; i < num_steps; ++i) {
        Step_sptr step(new Step(step_length));
        step->append(get_step("step", lattice_it, left, lattice_end,
                step_length), 1.0);
        get_steps().push_back(step);
    }
    if (lattice_it != lattice_end) {
        throw(std::runtime_error(
                "internal error: Independent_stepper did not make it to the end of the lattice\n"));
    }
    this->lattice_simulator.construct_sliced_chef_beamline(
            extract_slices(get_steps()));
}

Independent_stepper::~Independent_stepper()
{

}

//Independent_stepper_elements
Independent_stepper_elements::Independent_stepper_elements(
        Lattice_simulator const& lattice_simulator, int steps_per_element) :
    lattice_simulator(lattice_simulator)
{
    if (steps_per_element < 1) {
        throw std::runtime_error(
                "Independent_stepper_elements: steps_per_element must be >= 1");
    }
    for (Lattice_elements::iterator it =
            this->lattice_simulator.get_lattice_sptr()->get_elements().begin(); it
            != this->lattice_simulator.get_lattice_sptr()->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        if (length == 0.0) {
            Independent_operator_sptr
                    ind_op(
                            new Independent_operator(
                                    "step",
                                    this->lattice_simulator.get_operation_extractor_map_sptr()));
            Lattice_element_slice_sptr slice(new Lattice_element_slice(*(*it)));
            ind_op->append_slice(slice);
            Step_sptr step(new Step(0.0));
            step->append(ind_op, 1.0);
            get_steps().push_back(step);
        } else {
            double step_length = length / steps_per_element;
            for (int i = 0; i < steps_per_element; ++i) {
                Independent_operator_sptr
                        ind_op(
                                new Independent_operator(
                                        "step",
                                        this->lattice_simulator.get_operation_extractor_map_sptr()));
                double left = i * step_length;
                double right = (i + 1) * step_length;
                Lattice_element_slice_sptr slice(new Lattice_element_slice(
                        *(*it), left, right));
                ind_op->append_slice(slice);
                Step_sptr step(new Step(step_length));
                step->append(ind_op, 1.0);
                get_steps().push_back(step);
            }
        }
    }
    this->lattice_simulator.construct_sliced_chef_beamline(
            extract_slices(get_steps()));
}

Independent_stepper_elements::~Independent_stepper_elements()
{

}


//Split_operator_stepper

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
            // The rest of the element fits in the half step
            Lattice_element_slice_sptr slice(new Lattice_element_slice(
                    *(*lattice_it), left, right));
            retval->append_slice(slice);
            length += (right - left);
            ++lattice_it;
            left = 0.0;
            if (floating_point_equal(length, half_step_length, tolerance)) {
                if ((lattice_it == lattice_end) || ((*lattice_it)->get_length()
                        != 0.0)) {
                    complete = true;
                }
            } else {
                if (lattice_it == lattice_end) {
                    throw(std::runtime_error(
                            "get_half_step stepped beyond end of lattice"));
                }
            }
        } else {
            // Need to take a portion of the element...
            bool end_within_error = false;
            double old_right = right;
            right = half_step_length - length + left;
            if ((old_right - right) < retval->get_slices().size() * tolerance) {
                // ... unless we are within an accumulated tolerance of the end
                right = old_right;
                end_within_error = true;
            }
            Lattice_element_slice_sptr slice(new Lattice_element_slice(
                    *(*lattice_it), left, right));
            retval->append_slice(slice);
            length += (right - left);
            if (end_within_error) {
                ++lattice_it;
                left = 0.0;
            } else {
                left = right;
            }
            complete = true;
        }
    }
    return retval;
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
        Step_sptr step(new Step(step_length));
        step->append(get_half_step("first_half", lattice_it, left, lattice_end,
                half_step_length), 0.5);
        for (Collective_operators::const_iterator coll_op_it =
                collective_operators.begin(); coll_op_it
                != collective_operators.end(); ++coll_op_it) {
            step->append(*coll_op_it, 1.0);
        }
        step->append(get_half_step("second_half", lattice_it, left,
                lattice_end, half_step_length), 0.5);
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

void
Split_operator_stepper_elements::construct(
        Collective_operators const& collective_operators,
        int steps_per_element)
{
    if (steps_per_element < 1) {
        throw std::runtime_error(
                "Split_operator_stepper_elements: steps_per_element must be >= 1");
    }

    for (Lattice_elements::iterator it =
            this->lattice_simulator.get_lattice_sptr()->get_elements().begin(); it
            != this->lattice_simulator.get_lattice_sptr()->get_elements().end(); ++it) {
        double length = (*it)->get_length();

        //zero-length element
        if (length == 0.0) {
            Independent_operator_sptr ind_op(new Independent_operator("step", this->lattice_simulator.get_operation_extractor_map_sptr()));
            Lattice_element_slice_sptr slice(new Lattice_element_slice(*(*it)));
            ind_op->append_slice(slice);
            Step_sptr step(new Step(0.0));
            step->append(ind_op, 1.0);
            get_steps().push_back(step);
        //else
        } else {
            double step_length = length / steps_per_element;
            double half_step_length = 0.5 * step_length;

            for (int i = 0; i < steps_per_element; ++i) {
                double left = i * step_length;
                double middle = left + half_step_length;
                double right = (i + 1) * step_length;

                //1st Half
                Independent_operator_sptr
                        ind_op(
                                new Independent_operator(
                                        "step",
                                        this->lattice_simulator.get_operation_extractor_map_sptr()));
                Lattice_element_slice_sptr slice_1st_half(new Lattice_element_slice(
                        *(*it), left, middle));
                ind_op->append_slice(slice_1st_half);
                Step_sptr step(new Step(step_length));
                step->append(ind_op, 0.5);

                //Collective Effects
                for (Collective_operators::const_iterator coll_op_it =
                        collective_operators.begin(); coll_op_it
                        != collective_operators.end(); ++coll_op_it) {
                    step->append(*coll_op_it, 1.0);
                }

                //2nd Half
                Lattice_element_slice_sptr slice_2nd_half(new Lattice_element_slice(
                        *(*it), middle, right));
                //slice(new Lattice_element_slice(*(*it), middle, right));
                ind_op->append_slice(slice_2nd_half);
                step->append(ind_op, 0.5);

                get_steps().push_back(step);

            }
        }
     lattice_simulator.construct_sliced_chef_beamline(
            extract_slices(get_steps()));
    }
}

Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_simulator const& lattice_simulator,
        Collective_operator_sptr collective_operator, int steps_per_element) :
    lattice_simulator(lattice_simulator)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, steps_per_element);
}

Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_simulator const& lattice_simulator,
        Collective_operators const& collective_operators, int steps_per_element) :
    lattice_simulator(lattice_simulator)
{
    construct(collective_operators, steps_per_element);
}

Split_operator_stepper_elements::~Split_operator_stepper_elements()
{

}


