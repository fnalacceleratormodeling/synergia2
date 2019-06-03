#include "split_operator_stepper_elements.h"

std::vector<Step>
Split_operator_stepper_elements::apply_impl(Lattice const & lattice) const
{
    if (steps_per_element < 1) 
    {
        throw std::runtime_error(
                "Split_operator_stepper_elements: steps_per_element must be >= 1");
    }

    std::vector<Step> steps;

    for (auto const & ele : lattice.get_elements())
    {
        double length = ele.get_length();

        //zero-length element
        if (length == 0.0) 
        {
            steps.emplace_back(0.0);
            steps.back()
                .append_independent("step", 1.0)
                .append_slice(ele);
        } 
        else 
        {
            double step_length = length / steps_per_element;
            double half_step_length = 0.5 * step_length;

            for (int i = 0; i < steps_per_element; ++i) 
            {
                double left = i * step_length;
                double middle = left + half_step_length;
                double right = (i + 1) * step_length;

                steps.emplace_back(step_length);

                // 1st Half
                steps.back()
                    .append_independent("first_half", 0.5)
                    .append_slice(ele, left, middle);

                // Collective Effects
                // steps.append_collective();

                // 2nd Half
                steps.back()
                    .append_independent("second_half", 0.5)
                    .append_slice(ele, middle, right);
            }
        }
    }

    return steps;
}


#if 0
void
Split_operator_stepper_elements::construct(
        Collective_operators const& collective_operators, 
        int steps_per_element )
{
    if (steps_per_element < 1) 
    {
        throw std::runtime_error(
                "Split_operator_stepper_elements: steps_per_element must be >= 1");
    }

    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it) 
    {
        double length = (*it)->get_length();

        //zero-length element
        if (length == 0.0) {
            Independent_operator_sptr
                    ind_op(
                            new Independent_operator(
                                    "step",
                                    get_lattice_simulator().get_operation_extractor_map_sptr(),
                                    get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
            Lattice_element_slice_sptr slice(new Lattice_element_slice(*it));
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

                Step_sptr step(new Step(step_length));

                //1st Half
                Independent_operator_sptr
                        ind_op_first_half(
                                new Independent_operator(
                                        "first_half",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                Lattice_element_slice_sptr slice_1st_half(
                        new Lattice_element_slice(*it, left, middle));
                ind_op_first_half->append_slice(slice_1st_half);
                step->append(ind_op_first_half, 0.5);

                //Collective Effects
                for (Collective_operators::const_iterator coll_op_it =
                        collective_operators.begin(); coll_op_it
                        != collective_operators.end(); ++coll_op_it) {
                    Collective_operator_sptr copied_collective_operator_sptr(
                                                            (*coll_op_it)->clone());
                    step->append(copied_collective_operator_sptr, 1.0);
                }

                //2nd Half
                Independent_operator_sptr
                        ind_op_second_half(
                                new Independent_operator(
                                        "second_half",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                Lattice_element_slice_sptr slice_2nd_half(
                        new Lattice_element_slice(*it, middle, right));
                //slice(new Lattice_element_slice(*it, middle, right));
                ind_op_second_half->append_slice(slice_2nd_half);
                step->append(ind_op_second_half, 0.5);

                get_steps().push_back(step);

            }
        }
    }
    get_lattice_simulator().set_slices(extract_slices(get_steps()));
}

Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_sptr lattice_sptr, int map_order,
        Collective_operator_sptr collective_operator, int steps_per_element) 
    : Stepper(lattice_sptr, map_order)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, steps_per_element);
}
#endif

#if 0
Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_sptr lattice_sptr, int map_order,
        Collective_operators const& collective_operators, int steps_per_element) :
    Stepper(lattice_sptr, map_order)
{
    construct(collective_operators, steps_per_element);
}


#endif
