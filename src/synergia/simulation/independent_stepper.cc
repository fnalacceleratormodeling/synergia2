#include "independent_stepper.h"
#include "synergia/utils/string_utils.h"

void
Independent_stepper::construct(int num_steps)
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

Independent_stepper::Independent_stepper(
        Lattice_sptr lattice_sptr, int map_order, int num_steps) :
        Stepper(lattice_sptr, map_order)
{
    construct(num_steps);
}

Independent_stepper::Independent_stepper(
        Lattice_simulator const& lattice_simulator, int num_steps) :
                Stepper(lattice_simulator)
{
    construct(num_steps);
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
