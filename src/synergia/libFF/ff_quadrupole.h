#ifndef FF_QUADRUPOLE_H
#define FF_QUADRUPOLE_H

#include "ff_element.h"
#include "ff_drift.h"

class FF_quadrupole : public FF_element
{
private:
    static const int steps;
    static const int drifts_per_step;
    double get_reference_cdt(double length, double k,
                             Reference_particle & reference_particle);
public:
    FF_quadrupole();
    template <typename T>
    inline static void thin_quadrupole_unit(T const& x, T & xp,
                                            T const& y, T & yp, double const * kL);
#if 0
    template <typename T>
    inline static void thick_quadrupole_unit(T & x, T & xp,
                                             T & y, T & yp,
                                             T & cdt, T const& dpop,
                                             double reference_momentum,
                                             double m, double substep_reference_cdt,
                                             double step_length, double step_strength,
                                             int steps);
#endif

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_quadrupole();
};

typedef boost::shared_ptr<FF_quadrupole > FF_quadrupole_sptr;

template <typename T>
inline void FF_quadrupole::thin_quadrupole_unit(T const& x, T& xp,
                                                T const& y, T& yp, double const * kL) {
    xp += -kL[0] * x;
    yp += kL[0] * y;
}

#if 0
template <typename T>
inline void FF_quadrupole::thick_quadrupole_unit(T & x, T & xp,
                                         T & y, T & yp,
                                         T & cdt, T const& dpop,
                                         double reference_momentum,
                                         double m, double substep_reference_cdt,
                                         double step_length, double step_strength,
                                         int steps)
{
    // see yoshida4.py for formulas
    const double c1 = 0.675603595979828817023843904487;
    const double c4 = c1;
    const double c2 = -0.175603595979828817023843904487;
    const double c3 = c2;
    const double d1 = 1.35120719195965763404768780897;
    const double d3 = d1;
    const double d2 = -1.70241438391931526809537561795;

    for(int i = 0; i < steps; ++i) {
        FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                   m, substep_reference_cdt);
        thin_quadrupole_unit(x, xp, y, yp, d1 * step_strength);
        FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                   m, substep_reference_cdt);
        thin_quadrupole_unit(x, xp, y, yp, d2 * step_strength);
        FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                   m, substep_reference_cdt);
        thin_quadrupole_unit(x, xp, y, yp, d3 * step_strength);
        FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                   m, substep_reference_cdt);
    }
}
#endif

#endif // FF_QUADRUPOLE_H
