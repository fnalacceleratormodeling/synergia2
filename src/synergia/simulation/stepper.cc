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



Lattice_element_slices
Stepper::extract_slices(Steps const& steps)
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




