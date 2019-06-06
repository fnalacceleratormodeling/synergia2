#ifndef FF_DRIFT_H
#define FF_DRIFT_H

#include "synergia/libFF/ff_element.h"

class FF_drift : public FF_element
{

public:

    FF_drift() { }

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
};

#endif // FF_DRIFT_H
