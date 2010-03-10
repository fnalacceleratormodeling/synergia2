#ifndef STEPPER_H_
#define STEPPER_H_

#include <list>

#include "components/lattice/lattice.h"
#include "components/lattice/lattice_element_slice.h"
#include "components/lattice/chef_lattice.h"

class Operator
{
public:
    std::string name;
    Operator(std::string const& name)
    {
        this->name = name;
    }
    ;
    virtual void
    print() const
    {
        std::cout << "Operator\n";
    }
    virtual
    ~Operator()
    {
    }
    ;
};

typedef std::list<Operator > Operators;

class Collective_operator : public Operator
{
public:
    Collective_operator()
    {
    }
    ;
    virtual
    void
    print() const
    {
        std::cout << "Collective_operator\n";
    }
    ~Collective_operator()
    {
    }
    ;
};

typedef std::list<Collective_operator > Collective_operators;

typedef std::list<Lattice_element_slice > Lattice_element_slices;

class Independent_operator : public Operator
{
private:
    Lattice_element_slices slice_list;
public:
    Independent_operator();
    void
    append_slice(Lattice_element_slice const& slice);
    Lattice_element_slices const&
    get_slices() const;
    virtual void
    print() const;
    ~Independent_operator();
};

class Step
{
private:
    Operators operators;
public:
    Step();
    void
    append(Operator const& the_operator);
    void
    append(Operators const& the_operators);
    //    Operators const&
    //    get_operators() const;
    void
    print() const;
};

typedef std::list<Step > Steps;

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
    Independent_operator
    get_half_step(Lattice_element_list::iterator & lattice_it, double & left,
            Lattice_element_list::iterator const & lattice_end,
            const double half_step_length);
    void
    construct(Lattice & lattice, Chef_lattice & chef_lattice, int num_steps,
            Collective_operators const & collective_operators);
public:
    Split_operator_stepper(Lattice & lattice, Chef_lattice & chef_lattice,
            int num_steps, Collective_operator const & collective_operator);
    Split_operator_stepper(Lattice & lattice, Chef_lattice & chef_lattice,
            int num_steps, Collective_operators const & collective_operators);

};

#endif /* STEPPER_H_ */
