#ifndef FF_NONLINEARLENS_H
#define FF_NONLINEARLENS_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_element.h"

class FF_nllens : public FF_element
{
private:
    double get_reference_cdt(double icnll, double kick, Reference_particle & reference_particle);
public:
    FF_nllens();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_nllens();
};

typedef boost::shared_ptr<FF_nllens > FF_nllens_sptr;

#endif // FF_NONLINEARLENS_H
