#ifndef STEPPER_H_
#define STEPPER_H_

#include <list>
#include <boost/shared_ptr.hpp>

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
        std::cout << "Operator " << name << std::endl;
    }
    virtual
    ~Operator()
    {
    }
    ;
};

typedef boost::shared_ptr<Operator > Operator_sptr;
typedef std::list<Operator_sptr > Operators;

class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name) :
        Operator(name)
    {
    }
    ;
    virtual
    void
    print() const
    {
        std::cout << "Collective_operator: " << name << std::endl;
    }
    ~Collective_operator()
    {
    }
    ;
};

typedef boost::shared_ptr<Collective_operator > Collective_operator_sptr;
typedef std::list<Collective_operator_sptr > Collective_operators;

typedef boost::shared_ptr<Lattice_element_slice > Lattice_element_slice_sptr;
typedef std::list<Lattice_element_slice_sptr > Lattice_element_slices;

class Independent_operator : public Operator
{
private:
    Lattice_element_slices slices;
public:
    Independent_operator(std::string const& name);
    void
    append_slice(boost::shared_ptr<Lattice_element_slice > slice);
    Lattice_element_slices const&
    get_slices() const;
    virtual void
    print() const;
    ~Independent_operator();
};

typedef boost::shared_ptr<Independent_operator > Independent_operator_sptr;
typedef std::list<Independent_operator_sptr > Independent_operators;

class Step
{
private:
    Operators operators;
public:
    Step();
    void
    append(Operator_sptr operator_sptr);
    void
    append(Operators const& operators);
    //    Operators const&
    //    get_operators() const;
    void
    print(int index) const;
};

typedef boost::shared_ptr<Step > Step_sptr;
typedef std::list<Step_sptr > Steps;

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
            Lattice_element_list::iterator & lattice_it, double & left,
            Lattice_element_list::iterator const & lattice_end,
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
