#ifndef FF_ELENS_H
#define FF_ELENS_H

#include "ff_element.h"

class FF_elens : public FF_element
{
public:
    FF_elens() { };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_elens();
};

typedef boost::shared_ptr<FF_elens > FF_elens_sptr;

#endif // FF_HKICKER_H
