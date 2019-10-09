#include "split_operator_stepper.h"
#include "synergia/utils/string_utils.h"

namespace
{
    const std::string force_diag_attr = "force_diagnostics";
    const double fixed_step_tolerance = 1.0e-8;

    // Return an Independent_operator for a half step, starting at the
    // lattice_element given by lattice_it at position left. Both lattice_it
    // and left are updated by the function.

    template <typename ITER>
    std::vector<Lattice_element_slice> get_fixed_step_slices( 
            ITER & it, ITER const& end,
            double & left, double step_length, double & offset_fudge, 
            bool end_on_force_diagnostics )
    {
        std::vector<Lattice_element_slice> slices;

        double length = offset_fudge;
        bool complete = false;

        while (!complete) 
        {
            double right = it->get_length();

            if (length + (right - left) - fixed_step_tolerance <= step_length) 
            {
                // The rest of the element fits in the half step
                Lattice_element_slice slice(*it, left, right);
                slices.push_back(slice);
                length += (right - left);

                if (end_on_force_diagnostics 
                        && slice.has_right_edge() 
                        && it->has_string_attribute(force_diag_attr)
                        && !false_string(it->get_string_attribute(force_diag_attr)) )
                {
                    complete = true;
                }

                ++it;
                left = 0.0;

                if (std::abs(length - step_length) < fixed_step_tolerance) 
                { 
                    if ((it == end) || (it->get_length() != 0.0)) 
                        complete = true;
                } 
                else 
                {
                    if (it == end)
                        throw std::runtime_error("get_step stepped beyond end of lattice");
                }
            } 
            else 
            {
                // Need to take a portion of the element...
                bool end_within_error = false;
                double old_right = right;
                right = step_length - length + left;

                if (std::abs(old_right - right) < fixed_step_tolerance) 
                {
                    // ... unless we are within an accumulated tolerance of the end
                    right = old_right;
                    end_within_error = true;
                }

                Lattice_element_slice slice(*it, left, right);
                slices.push_back(slice);
                length += (right - left);

                if (end_within_error) 
                {
                    ++it;
                    left = 0.0;
                } 
                else 
                {
                    left = right;
                }

                complete = true;
            }
        }

        offset_fudge = length - step_length;
        return slices;
    }
}



std::vector<Step>
Split_operator_stepper::apply_impl(Lattice const & lattice) const
{
    std::vector<Step> steps;
    auto col_op_ptr = std::shared_ptr<Operator>(co_ops->create_operator());

    double step_length = lattice.get_length() / num_steps;
    double half_step_length = 0.5 * step_length;

    auto lattice_it  = lattice.get_elements().begin();
    auto lattice_end = lattice.get_elements().end();

    double left = 0.0;
    double offset_fudge = 0.0;

    for (int i = 0; i < num_steps; ++i) 
    {
        auto substep_lattice_it = lattice_it;

        double substep_left = left;
        double substep_offset_fudge = offset_fudge;

        // slices for first half step
        auto fst_half_slices = get_fixed_step_slices(lattice_it, lattice_end,
                left, half_step_length, offset_fudge, false);

        // slices for second half step
        auto snd_half_slices = get_fixed_step_slices(lattice_it, lattice_end,
                left, half_step_length, offset_fudge, false);

        // merge first half slices and second half slices
        std::vector<Lattice_element_slice> all_slices;
        all_slices.reserve(fst_half_slices.size() + snd_half_slices.size());
        all_slices.insert(all_slices.end(), fst_half_slices.begin(), fst_half_slices.end());
        all_slices.insert(all_slices.end(), snd_half_slices.begin(), snd_half_slices.end());

        double substep_length = 0.0;
        double all_substeps_length = 0.0;
        bool   found_force = false;

        for (auto it = all_slices.begin(); it != all_slices.end(); ++it) 
        {
            substep_length += it->get_right() - it->get_left();

            // jfa: I don't know of a simpler way to skip the last element of a
            //      container
            auto tmp_it = it; ++tmp_it;
            if (tmp_it != all_slices.end()) 
            {
                if (it->has_right_edge() 
                        && it->get_lattice_element().has_string_attribute(force_diag_attr)
                        && !false_string(it->get_lattice_element().get_string_attribute(force_diag_attr)) )
                {
                    found_force = true;
                    all_substeps_length += substep_length;

                    // create step
                    steps.emplace_back(substep_length);

                    if (substep_length == 0.0) 
                    {
                        auto sub_slices = get_fixed_step_slices(
                                substep_lattice_it, lattice_end,
                                substep_left, 0.0, substep_offset_fudge, true);

                        auto & op = steps.back().append_independent("zero_length", 1.0);
                        for (auto const& s : sub_slices) op.append_slice(s);
                    } 
                    else 
                    {
#if 0
                        double half_substep_length = 0.5 * substep_length;
                        Independent_operator_sptr subfirst_half_op_sptr(
                                Stepper::get_fixed_step("first_half",
                                        substep_lattice_it, substep_left,
                                        lattice_end, half_substep_length,
                                        substep_offset_fudge, false));
                        Independent_operator_sptr subsecond_half_op_sptr(
                                Stepper::get_fixed_step("second_half",
                                        substep_lattice_it, substep_left,
                                        lattice_end, half_substep_length,
                                        substep_offset_fudge, true));
                        substep_sptr->append(subfirst_half_op_sptr, 0.5);
                        for (Collective_operators::const_iterator coll_op_it =
                                collective_operators.begin();
                                coll_op_it != collective_operators.end();
                                ++coll_op_it) {
                            Collective_operator_sptr copied_collective_operator_sptr(
                                    (*coll_op_it)->clone());
                            substep_sptr->append(copied_collective_operator_sptr, 1.0);
                        }
                        substep_sptr->append(subsecond_half_op_sptr, 0.5);
#endif
                    }

                    //get_steps().push_back(substep_sptr);
                    substep_length = 0.0;
                }
            } 
            else 
            {
#if 0
                if (found_force) 
                {
                    double remain_length = step_length - all_substeps_length;
                    if (remain_length <= fixed_step_tolerance) {
                        remain_length = 0.0;
                    }
                    Step_sptr remain_sptr(new Step(remain_length));
                    if (remain_length == 0.0) {
                        Independent_operator_sptr zero_len_op_sptr(
                                Stepper::get_fixed_step("zero_length",
                                        substep_lattice_it, substep_left,
                                        lattice_end, 0.0, substep_offset_fudge,
                                        true));
                        remain_sptr->append(zero_len_op_sptr, 1.0);

                    } else {
                        double half_remain_length = 0.5 * remain_length;
                        Independent_operator_sptr remain_first_half_op_sptr(
                                Stepper::get_fixed_step("first_half",
                                        substep_lattice_it, substep_left,
                                        lattice_end, half_remain_length,
                                        substep_offset_fudge, false));
                        Independent_operator_sptr remain_second_half_op_sptr(
                                Stepper::get_fixed_step("second_half",
                                        substep_lattice_it, substep_left,
                                        lattice_end, half_remain_length,
                                        substep_offset_fudge, true));
                        remain_sptr->append(remain_first_half_op_sptr, 0.5);
                        for (Collective_operators::const_iterator coll_op_it =
                                collective_operators.begin();
                                coll_op_it != collective_operators.end();
                                ++coll_op_it) {
                            Collective_operator_sptr copied_collective_operator_sptr(
                                    (*coll_op_it)->clone());
                            remain_sptr->append(copied_collective_operator_sptr, 1.0);
                        }
                        remain_sptr->append(remain_second_half_op_sptr, 0.5);
                    }
                    get_steps().push_back(remain_sptr);
                    if (substep_lattice_it != lattice_it) {
                        throw(std::runtime_error(
                                "internal error: Split_operator_stepper created an inconsistent force_diagnostics step"));
                    }
                }
#endif
            }
        }

        if (!found_force) 
        {
            // new step
            steps.emplace_back(step_length);

            // 1st half
            auto & op1 = steps.back().append_independent("first_half", 0.5);
            for(auto const& s : fst_half_slices) op1.append_slice(s);

            // collective
            steps.back().append(col_op_ptr);

            // 2nd half
            auto & op2 = steps.back().append_independent("second_half", 0.5);
            for(auto const& s : snd_half_slices) op2.append_slice(s);
        }
    }

    if (lattice_it != lattice_end) 
    {
        throw(std::runtime_error(
                "internal error: split_operator_stepper did not make it to the end of the lattice\n"));
    }

    return steps;
}



#if 0
void
Split_operator_stepper::construct(
        Collective_operators const& collective_operators, int num_steps)
{
    double step_length =
            get_lattice_simulator().get_lattice_sptr()->get_length()
                    / num_steps;
    double half_step_length = 0.5 * step_length;
    Lattice_elements::iterator lattice_it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
    Lattice_elements::iterator lattice_end =
            get_lattice_simulator().get_lattice_sptr()->get_elements().end();
    double left = 0.0;
    double offset_fudge = 0.0;
    for (int i = 0; i < num_steps; ++i) {
        Lattice_elements::iterator substep_lattice_it = lattice_it;
        double substep_left = left;
        double substep_offset_fudge = offset_fudge;

        Independent_operator_sptr first_half_op_sptr(
                Stepper::get_fixed_step("first_half", lattice_it, left,
                        lattice_end, half_step_length, offset_fudge, false));
        Independent_operator_sptr second_half_op_sptr(
                Stepper::get_fixed_step("second_half", lattice_it, left,
                        lattice_end, half_step_length, offset_fudge, false));

        double substep_length = 0.0;
        double all_substeps_length = 0.0;
        bool found_force = false;
        Lattice_element_slices all_slices(first_half_op_sptr->get_slices());
        Lattice_element_slices tmp_copy(second_half_op_sptr->get_slices());
        all_slices.splice(all_slices.end(), tmp_copy);
        for (Lattice_element_slices::const_iterator it = all_slices.begin();
                it != all_slices.end(); ++it) {
            substep_length += (*it)->get_right() - (*it)->get_left();
            // jfa: I don't know of a simpler way to skip the last element of a
            //      container
            Lattice_element_slices::const_iterator tmp_it(it);
            ++tmp_it;
            if (tmp_it != all_slices.end()) {
                if ((*it)->has_right_edge()
                        && (*it)->get_lattice_element().has_string_attribute(
                                Stepper::force_diagnostics_attribute)) {
                    if (!false_string(
                            (*it)->get_lattice_element().get_string_attribute(
                                    Stepper::force_diagnostics_attribute))) {
                        found_force = true;
                        all_substeps_length += substep_length;
                        Step_sptr substep_sptr(new Step(substep_length));
                        if (substep_length == 0.0) {
                            Independent_operator_sptr zero_len_op_sptr(
                                    Stepper::get_fixed_step("zero_length",
                                            substep_lattice_it, substep_left,
                                            lattice_end, 0.0,
                                            substep_offset_fudge, true));
                            substep_sptr->append(zero_len_op_sptr, 1.0);
                        } else {
                            double half_substep_length = 0.5 * substep_length;
                            Independent_operator_sptr subfirst_half_op_sptr(
                                    Stepper::get_fixed_step("first_half",
                                            substep_lattice_it, substep_left,
                                            lattice_end, half_substep_length,
                                            substep_offset_fudge, false));
                            Independent_operator_sptr subsecond_half_op_sptr(
                                    Stepper::get_fixed_step("second_half",
                                            substep_lattice_it, substep_left,
                                            lattice_end, half_substep_length,
                                            substep_offset_fudge, true));
                            substep_sptr->append(subfirst_half_op_sptr, 0.5);
                            for (Collective_operators::const_iterator coll_op_it =
                                    collective_operators.begin();
                                    coll_op_it != collective_operators.end();
                                    ++coll_op_it) {
                                Collective_operator_sptr copied_collective_operator_sptr(
                                        (*coll_op_it)->clone());
                                substep_sptr->append(copied_collective_operator_sptr, 1.0);
                            }
                            substep_sptr->append(subsecond_half_op_sptr, 0.5);
                        }
                        get_steps().push_back(substep_sptr);
                        substep_length = 0.0;
                    }
                }
            } else {
                if (found_force) {
                    double remain_length = step_length - all_substeps_length;
                    if (remain_length <= fixed_step_tolerance) {
                        remain_length = 0.0;
                    }
                    Step_sptr remain_sptr(new Step(remain_length));
                    if (remain_length == 0.0) {
                        Independent_operator_sptr zero_len_op_sptr(
                                Stepper::get_fixed_step("zero_length",
                                        substep_lattice_it, substep_left,
                                        lattice_end, 0.0, substep_offset_fudge,
                                        true));
                        remain_sptr->append(zero_len_op_sptr, 1.0);

                    } else {
                        double half_remain_length = 0.5 * remain_length;
                        Independent_operator_sptr remain_first_half_op_sptr(
                                Stepper::get_fixed_step("first_half",
                                        substep_lattice_it, substep_left,
                                        lattice_end, half_remain_length,
                                        substep_offset_fudge, false));
                        Independent_operator_sptr remain_second_half_op_sptr(
                                Stepper::get_fixed_step("second_half",
                                        substep_lattice_it, substep_left,
                                        lattice_end, half_remain_length,
                                        substep_offset_fudge, true));
                        remain_sptr->append(remain_first_half_op_sptr, 0.5);
                        for (Collective_operators::const_iterator coll_op_it =
                                collective_operators.begin();
                                coll_op_it != collective_operators.end();
                                ++coll_op_it) {
                            Collective_operator_sptr copied_collective_operator_sptr(
                                    (*coll_op_it)->clone());
                            remain_sptr->append(copied_collective_operator_sptr, 1.0);
                        }
                        remain_sptr->append(remain_second_half_op_sptr, 0.5);
                    }
                    get_steps().push_back(remain_sptr);
                    if (substep_lattice_it != lattice_it) {
                        throw(std::runtime_error(
                                "internal error: Split_operator_stepper created an inconsistent force_diagnostics step"));
                    }
                }
            }
        }
        if (!found_force) {
            Step_sptr step(new Step(step_length));
            step->append(first_half_op_sptr, 0.5);
            for (Collective_operators::const_iterator coll_op_it =
                    collective_operators.begin();
                    coll_op_it != collective_operators.end(); ++coll_op_it) {
                Collective_operator_sptr copied_collective_operator_sptr(
                        (*coll_op_it)->clone());
                step->append(copied_collective_operator_sptr, 1.0);
            }
            step->append(second_half_op_sptr, 0.5);
            get_steps().push_back(step);
        }

    }
    if (lattice_it != lattice_end) {
        throw(std::runtime_error(
                "internal error: split_operator_stepper did not make it to the end of the lattice\n"));
    }
    get_lattice_simulator().set_slices(extract_slices(get_steps()));
}

Split_operator_stepper::Split_operator_stepper(
        Lattice_sptr lattice_sptr, int map_order,
        Collective_operator_sptr collective_operator, int num_steps) :
    Stepper(lattice_sptr, map_order)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, num_steps);
}

Split_operator_stepper::Split_operator_stepper(
        Lattice_sptr lattice_sptr, int map_order,
        Collective_operators const& collective_operators, int num_steps) :
    Stepper(lattice_sptr, map_order)
{
    construct(collective_operators, num_steps);
}

Split_operator_stepper::Split_operator_stepper(
        Lattice_simulator const& lattice_simulator,
        Collective_operator_sptr collective_operator, int num_steps) :
    Stepper(lattice_simulator)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, num_steps);
}

Split_operator_stepper::Split_operator_stepper(
        Lattice_simulator const& lattice_simulator,
        Collective_operators const& collective_operators, int num_steps) :
    Stepper(lattice_simulator)
{
    construct(collective_operators, num_steps);
}

Split_operator_stepper::Split_operator_stepper()
{

}

template<class Archive>
    void
    Split_operator_stepper::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
    }

template
void
Split_operator_stepper::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Split_operator_stepper::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Split_operator_stepper::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Split_operator_stepper::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Split_operator_stepper::~Split_operator_stepper()
{

}

BOOST_CLASS_EXPORT_IMPLEMENT(Split_operator_stepper)
#endif
