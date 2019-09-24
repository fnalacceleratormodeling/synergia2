#ifndef FF_MULTIPOLE_H
#define FF_MULTIPOLE_H

#include "ff_element.h"

class FF_multipole : public FF_element
{
public:
    FF_multipole();

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_MULTIPOLE_H
