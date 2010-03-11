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
    Chef_lattice * chef_lattice_ptr;

    Steps &
    get_steps();
    Lattice &
    get_lattice();
    Chef_lattice &
    get_chef_lattice();
    virtual void
    print() const;

    virtual
    ~Stepper()
    {
    }
    ;
};

//class Independent_stepper : public Stepper
//{
//
//};

class Split_operator_stepper : public Stepper
{
private:
    Independent_operator_sptr
    get_half_step(std::string const& name,
            Lattice_elements::iterator & lattice_it, double & left,
            Lattice_elements::iterator const & lattice_end,
            const double half_step_length);
    void
    construct(Lattice & lattice, Chef_lattice & chef_lattice, int num_steps,
            Collective_operators const & collective_operators);
public:
    Split_operator_stepper(Lattice & lattice, Chef_lattice & chef_lattice,
            int num_steps, Collective_operator_sptr collective_operator);
    Split_operator_stepper(Lattice & lattice, Chef_lattice & chef_lattice,
            int num_steps, Collective_operators const & collective_operators);

};

#endif /* STEPPER_H_ */
