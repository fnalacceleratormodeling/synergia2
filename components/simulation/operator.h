#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>
#include <list>
#include <boost/shared_ptr.hpp>

#include "components/bunch/bunch.h"
#include "components/lattice/lattice_element_slice.h"
#include "components/lattice/chef_lattice.h"
#include "components/simulation/independent_operation.h"
#include "components/simulation/independent_params.h"

class Operator
{
public:
    std::string name;
    Lattice_element_slices slices;

    Operator(std::string const& name);
    std::string const&
    get_name() const;
    virtual Lattice_element_slices &
    get_slices();
    virtual
    void
    apply(Bunch & bunch, Chef_lattice & chef_lattice);
    virtual void
    print() const;
    virtual
    ~Operator();
};

typedef boost::shared_ptr<Operator > Operator_sptr;
typedef std::list<Operator_sptr > Operators;

class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name);
    virtual
    void
    print() const;
    ~Collective_operator()

    ;
};

typedef boost::shared_ptr<Collective_operator > Collective_operator_sptr;
typedef std::list<Collective_operator_sptr > Collective_operators;

class Independent_operator : public Operator
{
private:
    Independent_params * params_ptr;
    Independent_operations operations;
    bool
    have_operations;
    void
    update_operations(Chef_lattice & chef_lattice);
    bool
    need_update();
public:
    Independent_operator(std::string const& name, Independent_params * params_ptr);
    void
    append_slice(boost::shared_ptr<Lattice_element_slice > slice);
    Lattice_element_slices const&
    get_slices() const;
    virtual void
    apply(Bunch & bunch, Chef_lattice & chef_lattice);
    virtual void
    print() const;
    ~Independent_operator();
};

typedef boost::shared_ptr<Independent_operator > Independent_operator_sptr;
typedef std::list<Independent_operator_sptr > Independent_operators;

#endif /* OPERATOR_H_ */
