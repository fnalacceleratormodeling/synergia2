#ifndef FF_MCMLENS_H
#define FF_MCMLENS_H

#include "ff_element.h"

class FF_mcmlens : public FF_element
{
public:
    FF_mcmlens() { };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_elens();
};

typedef boost::shared_ptr<FF_mcmlens > FF_mcmlens_sptr;

#endif // FF_MCMLENS_H
