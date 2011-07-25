#ifndef STEPPER_H_
#define STEPPER_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"

class Stepper
{
private:
    Steps steps;

public:
    Steps &
    get_steps();
    virtual void
    print() const;

    virtual
    ~Stepper();
};

typedef boost::shared_ptr<Stepper > Stepper_sptr;

/// The Independent_stepper class generates evenly-spaced Independent_operator
/// steps through a Lattice. No collective effects are included.
class Independent_stepper : public Stepper
{
private:
    Lattice_simulator lattice_simulator;
    Independent_operator_sptr
    get_step(std::string const& name, Lattice_elements::iterator & lattice_it,
            double & left, Lattice_elements::iterator const & lattice_end,
            const double half_step_length, double & offset_fudge);
public:
    /// Construct an Independent_stepper
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param num_steps the number of steps to take in the Lattice
    Independent_stepper(Lattice_simulator const& lattice_simulator,
            int num_steps);

    ~Independent_stepper();

};

/// The Independent_stepper_elements class generates a constant number of
/// Independent_operator steps per thick element. Thin elements are assigned
/// a single step each. No collective effects are included.
class Independent_stepper_elements : public Stepper
{
private:
    Lattice_simulator lattice_simulator;
public:
    /// Construct an Independent_stepper
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param steps_per_element the number of steps per thick element
    Independent_stepper_elements(Lattice_simulator const& lattice_simulator,
            int steps_per_element);

    ~Independent_stepper_elements();
};

typedef boost::shared_ptr<Independent_stepper_elements >
        Independent_stepper_elements_sptr;

/// The Split_operator_stepper class generates evenly-spaced split-operator
/// steps through a Lattice. One or more collective effects are included per
/// step.
class Split_operator_stepper : public Stepper
{
private:
    Lattice_simulator lattice_simulator;
    Independent_operator_sptr
    get_half_step(std::string const& name,
            Lattice_elements::iterator & lattice_it, double & left,
            Lattice_elements::iterator const & lattice_end,
            const double half_step_length, double & offset_fudge);
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
    ~Split_operator_stepper();
};

typedef boost::shared_ptr<Split_operator_stepper > Split_operator_stepper_sptr;

/// The Split_operator_stepper_elements class generates a constant number of
/// split-operator steps per thick element. Thin elements are assigned
/// a single step each. One or more collective effects are included per
/// step.
class Split_operator_stepper_elements : public Stepper
{
private:
    Lattice_simulator lattice_simulator;
    void
    construct(Collective_operators const & collective_operators,
            int steps_per_element);
public:
    /// Construct a Split_operator_stepper_elements with a single Collective_operator
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operator the Collective_operator to apply in each step
    /// @param steps_per_element the number of steps per thick element
            Split_operator_stepper_elements(
                    Lattice_simulator const& lattice_simulator,
                    Collective_operator_sptr collective_operator,
                    int steps_per_element);

    /// Construct a Split_operator_stepper_elements with multiple Collective_operators
    /// @param lattice_simulator the Lattice_simulator for the Lattice
    /// @param collective_operators the set of Collective_operators to apply
    ///        in each step
    /// @param steps_per_element the number of steps per thick element
    Split_operator_stepper_elements(Lattice_simulator const& lattice_simulator,
            Collective_operators const & collective_operators,
            int steps_per_element);
    ~Split_operator_stepper_elements();
};

/// Generate steps through lattice based on envelope shape.
/// Includes collective effects.
//class Split_operator_stepper_smart : public Stepper
//{
//
//};

#endif /* STEPPER_H_ */
