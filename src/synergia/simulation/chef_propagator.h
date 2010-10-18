#ifndef CHEF_PROPAGATOR_H_
#define CHEF_PROPAGATOR_H_

#include "components/bunch/bunch.h"
#include "components/lattice/chef_lattice.h"

class Chef_propagator
{
private:
    Chef_elements chef_elements;

public:
    Chef_propagator(Chef_elements const& chef_elements);
    void
    apply(Bunch & bunch);
};

#endif /* CHEF_PROPAGATOR_H_ */
