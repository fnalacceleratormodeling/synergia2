#include "stepper.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/string_utils.h"
#include "synergia/bunch/bunch.h"
#include <cmath>

const std::string Stepper::force_diagnostics_attribute("force_diagnostics");
const double Stepper::fixed_step_tolerance = 1.0e-8;

Stepper::Stepper(Lattice_simulator const& lattice_simulator) :
    lattice_simulator(lattice_simulator)
{

}

Stepper::Stepper()
{

}

Lattice_simulator &
Stepper::get_lattice_simulator()
{
    return lattice_simulator;
}

Steps &
Stepper::get_steps()
{
    return steps;
}

void
Stepper::force_update_operations_no_collective()
{
    int total_num = 1;
    double real_num = 1.0;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch bunch(lattice_simulator.get_lattice_sptr()->get_reference_particle(),
            total_num, real_num, commxx_sptr);
    int verbosity = 0;
    Logger logger(0);
    double time_step = 1.0;
    for (int i = 0; i < 6; ++i) {
        bunch.get_local_particles()[0][i] = 0.0;
    }
    Diagnosticss dummy_diagnosticss;
    for (Steps::iterator sit = steps.begin(); sit != steps.end(); ++sit) {
        for (Operators::iterator oit = (*sit)->get_operators().begin();
                oit != (*sit)->get_operators().end(); ++oit) {
            if ((*oit)->get_type() == Independent_operator::type_name) {
                (*oit)->apply(bunch, time_step, **sit, verbosity,
                        dummy_diagnosticss, logger);
            }
        }
    }
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

template<class Archive>
    void
    Stepper::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(steps);
        ar & BOOST_SERIALIZATION_NVP(lattice_simulator);
    }

template
void
Stepper::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Stepper::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Stepper::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Stepper::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Stepper::~Stepper()
{

}

// Return an Independent_operator for a half step, starting at the
// lattice_element given by lattice_it at position left. Both lattice_it
// and left are updated by the function.
Independent_operator_sptr
Stepper::get_fixed_step(std::string const& name,
        Lattice_elements::iterator & lattice_it, double & left,
        Lattice_elements::iterator const & lattice_end,
        const double step_length, double & offset_fudge,
        bool end_on_force_diagnostics)
{
    Independent_operator_sptr retval(
            new Independent_operator(name,
                    get_lattice_simulator().get_operation_extractor_map_sptr(),
                    get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
    double length = offset_fudge;
    bool complete = false;
    while (!complete) {
        double right = (*lattice_it)->get_length();
        if (floating_point_leq(length + (right - left), step_length,
                fixed_step_tolerance)) {
            // The rest of the element fits in the half step
            Lattice_element_slice_sptr slice(
                    new Lattice_element_slice(*lattice_it, left, right));
            retval->append_slice(slice);
            length += (right - left);
            if (end_on_force_diagnostics && slice->has_right_edge()
                    && (*lattice_it)->has_string_attribute(
                            force_diagnostics_attribute)) {
                if (!false_string(
                        (*lattice_it)->get_string_attribute(
                                force_diagnostics_attribute))) {
                    complete = true;
                }
            }
            ++lattice_it;
            left = 0.0;
            if (floating_point_equal(length, step_length,
                    fixed_step_tolerance)) {
                if ((lattice_it == lattice_end)
                        || ((*lattice_it)->get_length() != 0.0)) {
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
            if ((old_right - right) < fixed_step_tolerance) {
                // ... unless we are within an accumulated tolerance of the end
                right = old_right;
                end_within_error = true;
            }
            Lattice_element_slice_sptr slice(
                    new Lattice_element_slice(*lattice_it, left, right));
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
    offset_fudge = length - step_length;
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
                        element_slices(
                                boost::static_pointer_cast<Independent_operator >(
                                        *o_it)->get_slices());
                all_slices.splice(all_slices.end(), element_slices);
            }
        }
    }
    return all_slices;
}



Independent_stepper::Independent_stepper(
        Lattice_simulator const& lattice_simulator, int num_steps) :
        Stepper(lattice_simulator)
{
    double step_length =
            get_lattice_simulator().get_lattice_sptr()->get_length()
                    / num_steps;

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

        Independent_operator_sptr full_step_op_sptr(
                Stepper::get_fixed_step("step", lattice_it, left, lattice_end,
                        step_length, offset_fudge, false));

        double substep_length = 0.0;
        double all_substeps_length = 0.0;
        bool found_force = false;
        for (Lattice_element_slices::const_iterator it =
                full_step_op_sptr->get_slices().begin();
                it != full_step_op_sptr->get_slices().end(); ++it) {
            substep_length += (*it)->get_right() - (*it)->get_left();
            // jfa: I don't know of a simpler way to skip the last element of a
            //      container
            Lattice_element_slices::const_iterator tmp_it(it);
            ++tmp_it;
            if (tmp_it != full_step_op_sptr->get_slices().end()) {
                if ((*it)->has_right_edge()
                        && (*it)->get_lattice_element().has_string_attribute(
                                Stepper::force_diagnostics_attribute)) {
                    if (!false_string(
                            (*it)->get_lattice_element().get_string_attribute(
                                    Stepper::force_diagnostics_attribute))) {
                        found_force = true;
                        all_substeps_length += substep_length;
                        Independent_operator_sptr substep_op_sptr(
                                Stepper::get_fixed_step("step",
                                        substep_lattice_it, substep_left,
                                        lattice_end, substep_length,
                                        substep_offset_fudge, true));
                        Step_sptr substep_sptr(new Step(substep_length));
                        substep_sptr->append(substep_op_sptr, 1.0);
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
                    Independent_operator_sptr remain_op_sptr(
                            Stepper::get_fixed_step("step", substep_lattice_it,
                                    substep_left, lattice_end, remain_length,
                                    substep_offset_fudge, false));
                    Step_sptr remain_step_sptr(new Step(remain_length));
                    remain_step_sptr->append(remain_op_sptr, 1.0);
                    get_steps().push_back(remain_step_sptr);
                    if (substep_lattice_it != lattice_it) {
                        throw(std::runtime_error(
                                "internal error: Independent_stepper created an inconsistent force_diagnostics step"));
                    }
                }
            }
        }
        if (!found_force) {
            Step_sptr step_sptr(new Step(step_length));
            step_sptr->append(full_step_op_sptr, 1.0);
            get_steps().push_back(step_sptr);
        }
    }
    if (lattice_it != lattice_end) {
        throw(std::runtime_error(
                "internal error: Independent_stepper did not make it to the end of the lattice\n"));
    }
    get_lattice_simulator().set_slices(extract_slices(get_steps()));
}

Independent_stepper::Independent_stepper()
{

}

template<class Archive>
    void
    Independent_stepper::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
    }

template
void
Independent_stepper::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Independent_stepper::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Independent_stepper::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Independent_stepper::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Independent_stepper::~Independent_stepper()
{

}

BOOST_CLASS_EXPORT_IMPLEMENT(Independent_stepper)

//Independent_stepper_elements
Independent_stepper_elements::Independent_stepper_elements(
        Lattice_simulator const& lattice_simulator, int steps_per_element) :
    Stepper(lattice_simulator)
{
    if (steps_per_element < 1) {
        throw std::runtime_error(
                "Independent_stepper_elements: steps_per_element must be >= 1");
    }
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it) {
        double length = (*it)->get_length();
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
        } else {
            double step_length = length / steps_per_element;
            for (int i = 0; i < steps_per_element; ++i) {
                Independent_operator_sptr
                        ind_op(
                                new Independent_operator(
                                        "step",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                double left = i * step_length;
                double right = (i + 1) * step_length;
                Lattice_element_slice_sptr slice(
                        new Lattice_element_slice(*it, left, right));
                ind_op->append_slice(slice);
                Step_sptr step(new Step(step_length));
                step->append(ind_op, 1.0);
                get_steps().push_back(step);
            }
        }
    }
    get_lattice_simulator().set_slices(extract_slices(get_steps()));
}

Independent_stepper_elements::Independent_stepper_elements()
{

}


template<class Archive>
    void
    Independent_stepper_elements::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
    }

template
void
Independent_stepper_elements::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Independent_stepper_elements::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Independent_stepper_elements::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Independent_stepper_elements::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Independent_stepper_elements::~Independent_stepper_elements()
{

}

BOOST_CLASS_EXPORT_IMPLEMENT(Independent_stepper_elements)

//Split_operator_stepper




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
                                substep_sptr->append(*coll_op_it, 1.0);
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
                            remain_sptr->append(*coll_op_it, 1.0);
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
                step->append(*coll_op_it, 1.0);
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

void
Split_operator_stepper_elements::construct(
        Collective_operators const& collective_operators, int steps_per_element)
{
    if (steps_per_element < 1) {
        throw std::runtime_error(
                "Split_operator_stepper_elements: steps_per_element must be >= 1");
    }

    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it) {
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
                    step->append(*coll_op_it, 1.0);
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
        get_lattice_simulator().set_slices(extract_slices(get_steps()));
    }
}

Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_simulator const& lattice_simulator,
        Collective_operator_sptr collective_operator, int steps_per_element) :
    Stepper(lattice_simulator)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, steps_per_element);
}

Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_simulator const& lattice_simulator,
        Collective_operators const& collective_operators, int steps_per_element) :
    Stepper(lattice_simulator)
{
    construct(collective_operators, steps_per_element);
}

Split_operator_stepper_elements::Split_operator_stepper_elements()
{

}


template<class Archive>
    void
    Split_operator_stepper_elements::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Stepper);
    }

template
void
Split_operator_stepper_elements::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Split_operator_stepper_elements::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Split_operator_stepper_elements::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Split_operator_stepper_elements::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Split_operator_stepper_elements::~Split_operator_stepper_elements()
{

}

void
Split_operator_stepper_choice::construct_per_element_else()
{
    bool verbose=false;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t_turn, t_turn1;
    t_turn= MPI_Wtime();
    // check if all num_steps are >1
    for(List_choice_map::const_iterator lch_it=this->list_choice_map.begin(); lch_it!= this->list_choice_map.end();
                 ++lch_it){
            if (lch_it->second.num_steps < 1) {
                 std::cout<<"element name="<<  lch_it->first<<" num steps= "<<lch_it->second.num_steps<<std::endl;
                throw std::runtime_error(
                    "Split_operator_stepper_choice: all num_steps must be >= 1");
        }
    }

    if (this->num_steps_else<1) throw std::runtime_error(
                    "Split_operator_stepper_choice: num_steps_else must be >= 1");
    // find the length of the lattice without the chosen elements
    double length_else=0.;
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){

        List_choice_map::const_iterator lch_it= this->list_choice_map.find((*it)->get_name());
        if (lch_it == this->list_choice_map.end()) length_else +=(*it)->get_length();

    }

    double step_length_else = length_else/this->num_steps_else;

    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){
        double length = (*it)->get_length();
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
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" zero length step"<<std::endl;
         }
         else{
            List_choice_map::const_iterator lch_it= this->list_choice_map.find((*it)->get_name());
            if  (lch_it != this->list_choice_map.end()){
               double step_length= length/lch_it->second.num_steps;
               double half_step_length = 0.5 * step_length;
               int steps_per_element=lch_it->second.num_steps;
                for (int i = 0; i <steps_per_element;  ++i) {
                    double left = i * step_length;
                    double middle = left + half_step_length;
                    double right = (i + 1) * step_length;

                    Step_sptr step(new Step(step_length));
                    //1st Half
                    Independent_operator_sptr
                        ind_op_first_half(
                                new Independent_operator(
                                        "first_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_1st_half(
                            new Lattice_element_slice(*it, left, middle));
                    ind_op_first_half->append_slice(slice_1st_half);
                    step->append(ind_op_first_half, 0.5);


                     //Collective Effects
                    for (Collective_operators::const_iterator coll_op_it =
                        lch_it->second.collective_operators.begin(); coll_op_it
                        != lch_it->second.collective_operators.end();
                        ++coll_op_it) {
                    step->append(*coll_op_it, 1.0);
                    }

                        //2nd Half
                    Independent_operator_sptr
                        ind_op_second_half(
                                new Independent_operator(
                                        "second_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_2nd_half(
                        new Lattice_element_slice(*it, middle, right));
                    ind_op_second_half->append_slice(slice_2nd_half);
                    step->append(ind_op_second_half, 0.5);

                    get_steps().push_back(step);

                }
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" nr of collective kicks ="<<steps_per_element<<std::endl;
            }
            else{
                        int steps_per_element=int(ceil(length/step_length_else));
                        double step_length=length/steps_per_element;
                        double half_step_length=0.5*step_length;

                        for(int i = 0; i < steps_per_element; ++i){
                            double left = i * step_length;
                            double middle = left + half_step_length;
                            double right = (i + 1) * step_length;

                            Step_sptr step(new Step(step_length));
                          //1st Half
                            Independent_operator_sptr
                                ind_op_first_half(
                                    new Independent_operator(
                                            "first_half_else",
                                                get_lattice_simulator().get_operation_extractor_map_sptr(),
                                                    get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                            Lattice_element_slice_sptr slice_1st_half(
                                    new Lattice_element_slice(*it, left, middle));
                            ind_op_first_half->append_slice(slice_1st_half);
                            step->append(ind_op_first_half, 0.5);

                            //Collective Effects
                            for (Collective_operators::const_iterator coll_op_it =
                                this->list_choice_map["else"].collective_operators.begin(); coll_op_it
                                != this->list_choice_map["else"].collective_operators.end(); ++coll_op_it) {
                                step->append(*coll_op_it, 1.0);
                            }
                              //2nd Half
                            Independent_operator_sptr
                                ind_op_second_half(
                                    new Independent_operator(
                                            "second_half_else",
                                            get_lattice_simulator().get_operation_extractor_map_sptr(),
                                            get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                            Lattice_element_slice_sptr slice_2nd_half(
                                    new Lattice_element_slice(*it, middle, right));

                            ind_op_second_half->append_slice(slice_2nd_half);
                            step->append(ind_op_second_half, 0.5);

                            get_steps().push_back(step);

                   }
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" nr of \"else\" collective kicks ="<<steps_per_element<<std::endl;
            }
         }

    }

    get_lattice_simulator().set_slices(extract_slices(get_steps()));

    t_turn1= MPI_Wtime();
     if (rank==0){
          std::cout<<" stepper_choice per element:"<<std::endl;
          std::cout<<"                total number of steps=  "<<this->get_steps().size()<<std::endl;
          std::cout<<" time: stepper choice done in "<<t_turn1-t_turn<<std::endl;
     }

}

BOOST_CLASS_EXPORT_IMPLEMENT(Split_operator_stepper_elements);




void
Split_operator_stepper_choice::
make_stepper_else(Lattice_elements::iterator const& begin, Lattice_elements::iterator const& end,
    double const & length_between, double const & max_step_length)
{
    int num_steps=int(ceil(length_between/max_step_length));
    double step_length=length_between/num_steps;
    double  half_step_length=0.5*step_length;
    double left = 0.0;
    double offset_fudge = 0.0;
    Lattice_elements::iterator lattice_it =begin;
   // Lattice_elements::iterator lattice_end =end;
    for (int i = 0; i < num_steps; ++i) {
        Step_sptr step(new Step(step_length));
        step->append(
                Stepper::get_fixed_step("first_half_else", lattice_it, left, end,
                       half_step_length, offset_fudge, false), 0.5);



        for (Collective_operators::const_iterator coll_op_it =
                this->list_choice_map["else"].collective_operators.begin(); coll_op_it
                != this->list_choice_map["else"].collective_operators.end(); ++coll_op_it) {
            step->append(*coll_op_it, 1.0);
        }
        step->append(
                Stepper::get_fixed_step("second_half_else", lattice_it, left, end,
                        half_step_length, offset_fudge, false), 0.5);
        get_steps().push_back(step);
    }
    if (lattice_it != end) {
        throw(std::runtime_error(
                "internal error: split_operator_stepper_choice: make_stepper_else: did not make it to the end of the lattice\n"));
    }
}



void
Split_operator_stepper_choice::construct_split_else()
{
    bool verbose=false;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t_turn, t_turn1;
    t_turn= MPI_Wtime();
    // check if all num_steps are >1
    for(List_choice_map::const_iterator lch_it=this->list_choice_map.begin(); lch_it!= this->list_choice_map.end();
                 ++lch_it){
            if (lch_it->second.num_steps < 1) {
                 std::cout<<"element name="<<  lch_it->first<<" num steps= "<<lch_it->second.num_steps<<std::endl;
                throw std::runtime_error(
                    "Split_operator_stepper_choice: all num_steps must be >= 1");
        }
    }

    if (this->num_steps_else<1) throw std::runtime_error(
                    "Split_operator_stepper_choice: num_steps_else must be >= 1");

    // find the length of the lattice without the chosen elements
    double length_else=0.;
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){

        List_choice_map::const_iterator lch_it= this->list_choice_map.find((*it)->get_name());
        if (lch_it == this->list_choice_map.end()) length_else +=(*it)->get_length();

    }


    double step_length_else = length_else/this->num_steps_else;
    if (rank==0) std::cout<<" length_else="<< length_else<<" step_length_else="<< step_length_else<<std::endl;

    double  length_between=0.; // length between chosen elements
    Lattice_elements::iterator begin=get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
    Lattice_elements::iterator end=get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
    for (Lattice_elements::iterator it =
            get_lattice_simulator().get_lattice_sptr()->get_elements().begin(); it
            != get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++it){

        double length = (*it)->get_length();
        List_choice_map::const_iterator found_it= this->list_choice_map.find((*it)->get_name());
        if  (found_it != this->list_choice_map.end()){ // this is a chosen element
            // first make the step for "else" if length_else (i.e. the length between chosen elements) is nonzero
            if (begin != end) {
                make_stepper_else(begin, end, length_between, step_length_else);
            }
            if (length == 0.0) {
                Independent_operator_sptr
                        ind_op(
                                new Independent_operator(
                                        "step_zero_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                Lattice_element_slice_sptr slice(new Lattice_element_slice(*it));
                ind_op->append_slice(slice);
                Step_sptr step(new Step(0.0));
                step->append(ind_op, 1.0);
                get_steps().push_back(step);
                if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" zero length step"<<std::endl;
            }
            else{
               int steps_per_element=found_it->second.num_steps;
               double step_length= length/steps_per_element;
               double half_step_length = 0.5 * step_length;

                for (int i = 0; i <steps_per_element;  ++i) {
                    double left = i * step_length;
                    double middle = left + half_step_length;
                    double right = (i + 1) * step_length;

                    Step_sptr step(new Step(step_length));
                    //1st Half
                    Independent_operator_sptr
                        ind_op_first_half(
                                new Independent_operator(
                                        "first_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_1st_half(
                            new Lattice_element_slice(*it, left, middle));
                    ind_op_first_half->append_slice(slice_1st_half);
                    step->append(ind_op_first_half, 0.5);


                     //Collective Effects
                    for (Collective_operators::const_iterator coll_op_it =
                        found_it->second.collective_operators.begin(); coll_op_it
                        != found_it->second.collective_operators.end();
                        ++coll_op_it) {
                    step->append(*coll_op_it, 1.0);
                    }

                        //2nd Half
                    Independent_operator_sptr
                        ind_op_second_half(
                                new Independent_operator(
                                        "second_half_chosen",
                                        get_lattice_simulator().get_operation_extractor_map_sptr(),
                                        get_lattice_simulator().get_aperture_operation_extractor_map_sptr()));
                    Lattice_element_slice_sptr slice_2nd_half(
                        new Lattice_element_slice(*it, middle, right));
                    ind_op_second_half->append_slice(slice_2nd_half);
                    step->append(ind_op_second_half, 0.5);

                    get_steps().push_back(step);

                }
            if ((rank==0) && verbose) std::cout<<" element: "<<(*it)->get_name()<<" nr of collective kicks ="<<steps_per_element<<std::endl;
            }
        ++end;
        begin=end;
        length_between=0.;
        }
        else{
                 ++end;
                 length_between +=length;
        }
    }

    if (begin != end) {
                make_stepper_else(begin, end, length_between, step_length_else);
     }
     get_lattice_simulator().set_slices(extract_slices(get_steps()));


    t_turn1= MPI_Wtime();
     if (rank==0){
          std::cout<<" stepper_choice split_else:"<<std::endl;
          std::cout<<"                total number of steps=  "<<this->get_steps().size()<<std::endl;
          std::cout<<" time: stepper choice done in "<<t_turn1-t_turn<<std::endl;
     }

}



Split_operator_stepper_choice::Split_operator_stepper_choice(
                    Lattice_simulator const& lattice_simulator, List_choice_map const & list_choice_map, bool split_else):
Stepper(lattice_simulator),    list_choice_map(list_choice_map)
{

   try{
        (this->list_choice_map.find("else")!= this->list_choice_map.end()) ?
            this->num_steps_else=this->list_choice_map["else"].num_steps:
            throw std::runtime_error(
                    "Split_operator_stepper_choice: if there is no  keyword \"else\" in the list_choice_map"
                     " you should use the constructor which provides \"num_steps_else\"");

        split_else ? construct_split_else():  construct_per_element_else();
    }
    catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 222);
    }
}

Split_operator_stepper_choice::Split_operator_stepper_choice(int num_steps_else,
                    Lattice_simulator const& lattice_simulator, List_choice_map const & list_choice_map,  bool split_else):
Stepper(lattice_simulator),    list_choice_map(list_choice_map), num_steps_else(num_steps_else)
{

    try{
        split_else ? construct_split_else():  construct_per_element_else();
    }
    catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 222);
    }

}

Split_operator_stepper_choice::~Split_operator_stepper_choice()
{
}
