#include "split_operator_stepper_elements.h"

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
        Collective_operator_sptr collective_operator, int steps_per_element) :
    Stepper(lattice_sptr, map_order)
{
    Collective_operators collective_operators;
    collective_operators.push_back(collective_operator);
    construct(collective_operators, steps_per_element);
}

Split_operator_stepper_elements::Split_operator_stepper_elements(
        Lattice_sptr lattice_sptr, int map_order,
        Collective_operators const& collective_operators, int steps_per_element) :
    Stepper(lattice_sptr, map_order)
{
    construct(collective_operators, steps_per_element);
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
BOOST_CLASS_EXPORT_IMPLEMENT(Split_operator_stepper_elements);


