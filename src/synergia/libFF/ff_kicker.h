#ifndef FF_KICKER_H
#define FF_KICKER_H

#include "ff_element.h"

class FF_kicker : public FF_element
{
private:
    double get_reference_cdt(double length, double hk, double vk, 
            Reference_particle & reference_particle);

public:
    FF_kicker() { };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_kicker();
};

typedef boost::shared_ptr<FF_kicker > FF_kicker_sptr;

#endif // FF_KICKER_H
