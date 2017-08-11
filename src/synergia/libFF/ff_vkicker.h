#ifndef FF_VKICKER_H
#define FF_VKICKER_H

#include "ff_element.h"

class FF_vkicker : public FF_element
{
private:
    double get_reference_cdt(double length, double k, Reference_particle & reference_particle);

public:
    FF_vkicker() { };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_vkicker();
};

typedef boost::shared_ptr<FF_vkicker > FF_vkicker_sptr;

#endif // FF_VKICKER_H
