#ifndef SPLIT_OPERATOR_STEPPER_ELEMENTS_H_
#define SPLIT_OPERATOR_STEPPER_ELEMENTS_H_
#include "synergia/simulation/stepper.h"

/// The Split_operator_stepper_elements class generates a constant number of
/// split-operator steps per thick element. Thin elements are assigned
/// a single step each. One or more collective effects are included per
/// step.
class Split_operator_stepper_elements : public Stepper
{
private:
    void
    construct(Collective_operators const & collective_operators,
            int steps_per_element);
public:
    /// Construct a Split_operator_stepper_elements with a single Collective_operator
    /// @param lattice_sptr the Lattice
    /// @param map_order order for Chef_map operations
    /// @param collective_operator the Collective_operator to apply in each step
    /// @param steps_per_element the number of steps per thick element
    Split_operator_stepper_elements(Lattice_sptr lattice_sptr, int map_order,
            Collective_operator_sptr collective_operator,
            int steps_per_element);

    /// Construct a Split_operator_stepper_elements with multiple Collective_operators
    /// @param lattice_sptr the Lattice
    /// @param map_order order for Chef_map operations
    /// @param collective_operators the set of Collective_operators to apply
    ///        in each step
    /// @param steps_per_element the number of steps per thick element
    Split_operator_stepper_elements(Lattice_sptr lattice_sptr, int map_order,
            Collective_operators const & collective_operators,
            int steps_per_element);

    /// Deprecated. Construct a Split_operator_stepper_elements with a single Collective_operator
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operator the Collective_operator to apply in each step
    /// @param steps_per_element the number of steps per thick element
    Split_operator_stepper_elements(Lattice_simulator const& lattice_simulator,
            Collective_operator_sptr collective_operator,
            int steps_per_element);

    /// Deprecated. Construct a Split_operator_stepper_elements with multiple Collective_operators
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operators the set of Collective_operators to apply
    ///        in each step
    /// @param steps_per_element the number of steps per thick element
    Split_operator_stepper_elements(Lattice_simulator const& lattice_simulator,
            Collective_operators const & collective_operators,
            int steps_per_element);

    /// Default constructor for serialization use only
    Split_operator_stepper_elements();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Split_operator_stepper_elements();
};
BOOST_CLASS_EXPORT_KEY(Split_operator_stepper_elements);
typedef boost::shared_ptr<Split_operator_stepper_elements >
        Split_operator_stepper_elements_sptr;  // syndoc:include

#endif /* SPLIT_OPERATOR_STEPPER_ELEMENTS_H_ */
