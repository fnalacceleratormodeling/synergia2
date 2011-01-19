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

/// Generate evenly-spaced steps through lattice with collective effects.
//class Independent_stepper : public Stepper
//{
//private:
//    Lattice_simulator lattice_simulator;
//    void
//    construct(Collective_operators const & collective_operators, int num_steps);
//public:
//    Independent_stepper(Lattice_simulator const& lattice_simulator,
//            Collective_operator_sptr collective_operator, int num_steps);
//    Independent_stepper(Lattice_simulator const& lattice_simulator,
//            Collective_operators const & collective_operators, int num_steps);
//    ~Independent_stepper();
//
//};

/// Generate per-element steps through lattice without collective effects.
class Independent_stepper_elements : public Stepper
{
private:
    Lattice_simulator lattice_simulator;
public:
    Independent_stepper_elements(Lattice_simulator const& lattice_simulator,
            int steps_per_element);
    ~Independent_stepper_elements();
};

typedef boost::shared_ptr<Independent_stepper_elements >
        Independent_stepper_elements_sptr;

/// Generate evenly-spaced steps through lattice with collective effects.
class Split_operator_stepper : public Stepper
{
private:
    Lattice_simulator lattice_simulator;
    Independent_operator_sptr
    get_half_step(std::string const& name,
            Lattice_elements::iterator & lattice_it, double & left,
            Lattice_elements::iterator const & lattice_end,
            const double half_step_length);
    void
    construct(Collective_operators const & collective_operators, int num_steps);
public:
    Split_operator_stepper(Lattice_simulator const& lattice_simulator,
            Collective_operator_sptr collective_operator, int num_steps);
    Split_operator_stepper(Lattice_simulator const& lattice_simulator,
            Collective_operators const & collective_operators, int num_steps);
    ~Split_operator_stepper();
};

typedef boost::shared_ptr<Split_operator_stepper > Split_operator_stepper_sptr;

/// Generate per-element steps through lattice with collective effects.
//class Split_operator_stepper_elements : public Stepper
//{
//
//};

/// Generate steps through lattice based on envelope shape.
/// Includes collective effects.
//class Split_operator_stepper_smart : public Stepper
//{
//
//};

#endif /* STEPPER_H_ */
