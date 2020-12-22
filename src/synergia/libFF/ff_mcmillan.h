#ifndef FF_MCMILLAN_H
#define FF_MCMILLAN_H

#include "ff_element.h"

class FF_mcmillan : public FF_element
{
private:
    double get_reference_cdt(double length, double j0, double beta_e, double radius, Reference_particle& rp, bool simple);

public:
    FF_mcmillan() { };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_mcmillan();
};

typedef boost::shared_ptr<FF_mcmillan > FF_mcmillan_sptr;

#endif // FF_MCMILLAN_H
