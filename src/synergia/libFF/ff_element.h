#ifndef FF_ELEMENT_H
#define FF_ELEMENT_H

#include <beamline/JetParticle.h>
#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element_slice.h"

class FF_element
{
public:
    FF_element();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle) = 0;
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch) = 0;
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_element();
};

typedef boost::shared_ptr<FF_element > FF_element_sptr; // syndoc:include

#endif // FF_ELEMENT_H
