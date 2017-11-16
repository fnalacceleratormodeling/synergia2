#ifndef FF_DIPEDGE_H
#define FF_DIPEDGE_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_element.h"

class FF_dipedge : public FF_element
{
private:
    double get_reference_cdt(double re_2_1, double re_4_3, double * te, Reference_particle & reference_particle);
public:
    FF_dipedge();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_dipedge();
};

typedef boost::shared_ptr<FF_dipedge > FF_dipedge_sptr;

#endif // FF_DIPEDGE_H
