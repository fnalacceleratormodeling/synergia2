#ifndef FF_SEXTUPOLE_H
#define FF_SEXTUPOLE_H

#include "ff_element.h"

class FF_sextupole : public FF_element
{
public:
    FF_sextupole() { /*steps = 1;*/ };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_SEXTUPOLE_H
