#ifndef STEPPER_H_
#define STEPPER_H_

#include <list>
#include <boost/shared_ptr.hpp>

#include "components/lattice/lattice.h"
#include "components/lattice/chef_lattice.h"

#include "components/simulation/operator.h"
#include "components/simulation/step.h"

class Stepper
{
public:
    Steps steps;
    Lattice * lattice_ptr;

    Steps &
    get_steps();
    Lattice &
    get_lattice();
    virtual void
    print() const;

    virtual
    ~Stepper();
};

/// Generate evenly-spaced steps through lattice without collective effects.
//class Independent_stepper : public Stepper
//{
//
//};

/// Generate per-element steps through lattice without collective effects.
//class Independent_stepper_elements: public Stepper
//{
//
//};

/// Generate evenly-spaced steps through lattice with collective effects.
class Split_operator_stepper : public Stepper
{
private:
    Independent_operator_sptr
    get_half_step(std::string const& name,
            Lattice_elements::iterator & lattice_it, double & left,
            Lattice_elements::iterator const & lattice_end,
            const double half_step_length);
    void
    construct(Lattice & lattice, int num_steps,
            Collective_operators const & collective_operators);
public:
    Split_operator_stepper(Lattice & lattice, int num_steps,
            Collective_operator_sptr collective_operator);
    Split_operator_stepper(Lattice & lattice, int num_steps,
            Collective_operators const & collective_operators);

};

/// Generate per-element steps through lattice with collective effects.
//class Split_operator_stepper_elements : public Stepper
//{
//
//};

/// Generate steps through lattice based on envelope shape.
/// Includes collective effects.
//class Split_operator_stepper_intelligent : public Stepper
//{
//
//};

#endif /* STEPPER_H_ */
