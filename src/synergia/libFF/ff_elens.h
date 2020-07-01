#ifndef FF_ELENS_H
#define FF_ELENS_H

#include "synergia/libFF/ff_element.h"

class FF_elens : public FF_element
{
public:
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_HKICKER_H
