#ifndef FF_OCTUPOLE_H
#define FF_OCTUPOLE_H

#include "ff_element.h"

class FF_octupole : public FF_element
{
public:
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_SEXTUPOLE_H
