#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>
#include <list>
#include <boost/shared_ptr.hpp>

#include "components/bunch/bunch.h"
#include "components/lattice/lattice_element_slice.h"
#include "components/lattice/chef_lattice.h"
#include "components/simulation/independent_operation.h"
#include "components/simulation/operation_extractor.h"

class Operator;
typedef boost::shared_ptr<Operator > Operator_sptr;
typedef std::list<Operator_sptr > Operators;

class Operator
{
private:
    std::string name, type;
public:
    Operator(std::string const& name, std::string const& type);
    std::string const&
    get_name() const;
    std::string const&
    get_type() const;
    virtual
    void
    apply(Bunch & bunch, Operators & step_operators) = 0;
    virtual void
    print() const;
    virtual
    ~Operator();
};


class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name);
    virtual
    void
    apply(Bunch & bunch, Operators & step_operators);
    ~Collective_operator();
};

typedef boost::shared_ptr<Collective_operator > Collective_operator_sptr;
typedef std::list<Collective_operator_sptr > Collective_operators;

class Independent_operator : public Operator
{
private:
    Lattice_element_slices slices;
    Independent_operations operations;
    Operation_extractor_map_sptr operation_extractor_map_sptr;
    bool have_operations;
    void
    update_operations(Reference_particle const& reference_particle);
    bool
    need_update();
public:
    Independent_operator(std::string const& name,
            Operation_extractor_map_sptr const& operation_extractor_map_sptr);
    void
    append_slice(Lattice_element_slice_sptr slice_sptr);
    Lattice_element_slices const&
    get_slices() const;
    virtual void
    apply(Bunch & bunch, Operators & step_operators);
    virtual void
    print() const;
    ~Independent_operator();
};

typedef boost::shared_ptr<Independent_operator > Independent_operator_sptr;
typedef std::list<Independent_operator_sptr > Independent_operators;

#endif /* OPERATOR_H_ */
