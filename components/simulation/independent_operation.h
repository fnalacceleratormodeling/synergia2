#ifndef INDEPENDENT_OPERATION_H_
#define INDEPENDENT_OPERATION_H_

#include "components/simulation/fast_mapping.h"
#include "components/simulation/chef_propagator.h"
#include "components/lattice/chef_lattice.h"
#include "boost/shared_ptr.hpp"
#include <list>
#include <string>
#include <map>

class Independent_operation
{
private:
    std::string type;
public:
    Independent_operation(std::string const& type);
    std::string
    get_type() const;
    virtual void
    apply(Bunch & bunch) = 0;
    virtual
    ~Independent_operation();
};

typedef boost::shared_ptr<Independent_operation > Independent_operation_sptr;
typedef std::list<Independent_operation_sptr > Independent_operations;

class Fast_mapping_operation : public Independent_operation
{
private:
    Fast_mapping mapping;

public:
    Fast_mapping_operation(Fast_mapping const& mapping);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Fast_mapping_operation();
};

typedef boost::shared_ptr<Fast_mapping_operation > Fast_mapping_operation_sptr;

class Chef_propagate_operation : public Independent_operation
{
private:
    Chef_propagator chef_propagator;

public:
    Chef_propagate_operation(Chef_elements const& elements);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Chef_propagate_operation();
};

typedef boost::shared_ptr<Chef_propagate_operation >
        Chef_propagate_operation_sptr;

#endif /* INDEPENDENT_OPERATION_H_ */
