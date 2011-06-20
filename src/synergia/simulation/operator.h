#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>
#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/simulation/independent_operation.h"
#include "synergia/simulation/operation_extractor.h"

class Step;

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
    virtual void
    apply(Bunch & bunch, double time_step, Step & step) = 0;
    virtual void
    print() const;
    ~Operator();
};
typedef boost::shared_ptr<Operator > Operator_sptr;
typedef std::list<Operator_sptr > Operators;

class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step) = 0;
    ~Collective_operator();
};

typedef boost::shared_ptr<Collective_operator > Collective_operator_sptr;
typedef std::list<Collective_operator_sptr > Collective_operators;

class Dummy_collective_operator : public Collective_operator
{
public:
    Dummy_collective_operator(std::string const& name);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step);
    ~Dummy_collective_operator();
};

typedef boost::shared_ptr<Dummy_collective_operator >
        Dummy_collective_operator_sptr;
typedef std::list<Dummy_collective_operator_sptr > Dummy_collective_operators;

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
            Operation_extractor_map_sptr operation_extractor_map_sptr);
    void
    append_slice(Lattice_element_slice_sptr slice_sptr);
    Lattice_element_slices const&
    get_slices() const;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step);
    virtual void
    print() const;
    ~Independent_operator();
};

typedef boost::shared_ptr<Independent_operator > Independent_operator_sptr;
typedef std::list<Independent_operator_sptr > Independent_operators;

#endif /* OPERATOR_H_ */
