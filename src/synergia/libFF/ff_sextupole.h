#ifndef FF_SEXTUPOLE_H
#define FF_SEXTUPOLE_H

#include "ff_element.h"
#include "ff_drift.h"

class FF_sextupole : public FF_element
{
private:
    static const int drifts_per_step;
    double get_reference_cdt(double length, double * k,
                             Reference_particle & reference_particle);
public:
    FF_sextupole();
    template <typename T>
    inline static void thin_sextupole_unit(T const& x, T & xp,
                                            T const& y, T & yp, double const * kL);

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_sextupole();
};

typedef boost::shared_ptr<FF_sextupole> FF_sextupole_sptr;

template <typename T>
inline void FF_sextupole::thin_sextupole_unit(T const& x, T& xp,
                                              T const& y, T& yp, double const * kL) {
    xp += -0.5 * kL[0] * (x * x - y * y) - kL[1] * x * y;
    yp += kL[0] * x * y - 0.5 * kL[1] * (x * x - y * y);
}

#endif // FF_SEXTUPOLE_H
