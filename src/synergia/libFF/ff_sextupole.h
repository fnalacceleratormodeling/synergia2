#ifndef FF_SEXTUPOLE_H
#define FF_SEXTUPOLE_H

#include "ff_element.h"

class FF_sextupole : public FF_element
{
private:
    double get_reference_cdt(double length, double * k,
                             Reference_particle & reference_particle);
public:
    FF_sextupole() { steps = 1; };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_sextupole();
};

typedef boost::shared_ptr<FF_sextupole> FF_sextupole_sptr;
#endif // FF_SEXTUPOLE_H
