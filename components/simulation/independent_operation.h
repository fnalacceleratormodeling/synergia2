#ifndef INDEPENDENT_OPERATION_H_
#define INDEPENDENT_OPERATION_H_

#include "components/simulation/fast_mapping.h"
#include "components/lattice/chef_lattice.h"
#include "boost/shared_ptr.hpp"
#include <list>
#include <string>
#include <map>

class Independent_operation
{
public:
    virtual void
    apply(Bunch & bunch) = 0;
    virtual
    ~Independent_operation()
    {
    }
    ;
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

class Chef_propagate_operation : public Independent_operation
{
private:
    Chef_elements chef_elements;

public:
    Chef_propagate_operation(Chef_elements const& elements);
    virtual void
    apply(Bunch & bunch);
    virtual
    ~Chef_propagate_operation();
};

#endif /* INDEPENDENT_OPERATION_H_ */
