#include "independent_stepper_elements.h"


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
