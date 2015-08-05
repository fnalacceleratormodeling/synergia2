#ifndef FF_HKICKER_H
#define FF_HKICKER_H

#include "ff_element.h"

class FF_hkicker : public FF_element
{
private:
    double get_reference_cdt(double length, double k, Reference_particle & reference_particle);

public:
    FF_hkicker() { };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_hkicker();
};

typedef boost::shared_ptr<FF_hkicker > FF_hkicker_sptr;

#endif // FF_HKICKER_H
