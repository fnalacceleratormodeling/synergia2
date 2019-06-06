#ifndef FF_QUADRUPOLE_H
#define FF_QUADRUPOLE_H

#include "ff_element.h"

class FF_quadrupole : public FF_element
{
public:
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_QUADRUPOLE_H
