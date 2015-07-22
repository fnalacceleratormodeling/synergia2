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
public:
    FF_rbend();

    template <typename T>
    inline static void thin_rbend_unit(T const & x, T & xp,
                                       T const & y, T & yp, double const * kL);

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_rbend();
};

template <typename T>
inline void FF_rbend::thin_rbend_unit(T const& x, T& xp,
                                      T const& y, T& yp, double const * kL) 
{
    FF_algorithm::thin_dipole_unit(x, xp, y, yp, kL);
    FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kL + 2);
    FF_algorithm::thin_sextupole_unit(x, xp, y, yp, kL + 4);
}


typedef boost::shared_ptr<FF_rbend > FF_rbend_sptr;

#endif // FF_RBEND_H
