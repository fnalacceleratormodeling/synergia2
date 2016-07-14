#ifndef FF_OCTUPOLE_H
#define FF_OCTUPOLE_H

#include "ff_element.h"

class FF_octupole : public FF_element
{
private:
    static const int drifts_per_step;
    double get_reference_cdt(double length, double * k,
                             Reference_particle & reference_particle);
public:
    FF_octupole();
    template <typename T>
    inline static void thin_octupole_unit(T const& x, T & xp,
                                            T const& y, T & yp, double const * kL);

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_octupole();
};

typedef boost::shared_ptr<FF_octupole> FF_octupole_sptr;
#endif // FF_SEXTUPOLE_H
