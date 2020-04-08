#ifndef FF_HKICKER_H
#define FF_HKICKER_H

#include "ff_element.h"

class FF_hkicker : public FF_element
{
public:
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_HKICKER_H
