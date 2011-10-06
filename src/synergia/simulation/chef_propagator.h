#ifndef CHEF_PROPAGATOR_H_
#define CHEF_PROPAGATOR_H_

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/lattice/chef_lattice_section.h"

class Chef_propagator
{
private:
    Chef_lattice_section_sptr chef_lattice_section_sptr;

public:
    Chef_propagator(Chef_lattice_section_sptr chef_lattice_section_sptr);
    void
    apply(Bunch & bunch);
};

#endif /* CHEF_PROPAGATOR_H_ */
