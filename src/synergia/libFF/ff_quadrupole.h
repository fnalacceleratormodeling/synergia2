#ifndef FF_QUADRUPOLE_H
#define FF_QUADRUPOLE_H

#include "ff_element.h"

class FF_quadrupole : public FF_element
{
private:
    double get_reference_cdt(double length, double * k,
                             Reference_particle & reference_particle);
public:
    FF_quadrupole() { /* order=4; steps=4; */ };

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_quadrupole();
};

typedef boost::shared_ptr<FF_quadrupole > FF_quadrupole_sptr;

#endif // FF_QUADRUPOLE_H
