#ifndef SPLIT_OPERATOR_STEPPER_H_
#define SPLIT_OPERATOR_STEPPER_H_
#include "synergia/simulation/stepper.h"

/// The Split_operator_stepper class generates evenly-spaced split-operator
/// steps through a Lattice. One or more collective effects are included per
/// step.
class Split_operator_stepper : public Stepper
{
    void
    construct(Collective_operators const & collective_operators, int num_steps);
public:
    /// Construct a Split_operator_stepper with a single Collective_operator
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operator the Collective_operator to apply in each step
    /// @param num_steps the number of steps to take in the Lattice
    Split_operator_stepper(Lattice_simulator const& lattice_simulator,
            Collective_operator_sptr collective_operator, int num_steps);

    /// Construct a Split_operator_stepper with multiple Collective_operators
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operators the set of Collective_operators to apply in each step
    /// @param num_steps the number of steps to take in the Lattice
    Split_operator_stepper(Lattice_simulator const& lattice_simulator,
            Collective_operators const & collective_operators, int num_steps);

    /// Default constructor for serialization use only
    Split_operator_stepper();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Split_operator_stepper();
};
BOOST_CLASS_EXPORT_KEY(Split_operator_stepper);
typedef boost::shared_ptr<Split_operator_stepper > Split_operator_stepper_sptr;


#endif /* SPLIT_OPERATOR_STEPPER_H_ */
