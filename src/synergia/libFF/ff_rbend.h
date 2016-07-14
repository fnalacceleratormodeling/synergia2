#ifndef FF_RBEND_H
#define FF_RBEND_H

#include "ff_element.h"
#include "ff_algorithm.h"
#include "synergia/utils/invsqrt.h"

class FF_rbend : public FF_element
{
private:

    double get_reference_cdt(double length, double * k,
                                   Reference_particle &reference_particle);

    double get_reference_cdt(double length, double strength, double angle,
                                   std::complex<double> const & phase,
                                   std::complex<double> const & term,
                                   Reference_particle &reference_particle);

public:
    FF_rbend();

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_rbend();
};

typedef boost::shared_ptr<FF_rbend > FF_rbend_sptr;

#endif // FF_RBEND_H
