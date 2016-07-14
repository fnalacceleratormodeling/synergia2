#ifndef FF_MULTIPOLE_H
#define FF_MULTIPOLE_H

#include "ff_element.h"

class FF_multipole : public FF_element
{
public:
    FF_multipole();

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_multipole();
};

typedef boost::shared_ptr<FF_multipole > FF_multipole_sptr;

#endif // FF_MULTIPOLE_H
