#ifndef FF_DRIFT_H
#define FF_DRIFT_H

#include "synergia/libFF/ff_element.h"

class FF_drift : public FF_element
{
private:
    double get_reference_cdt(double length, Reference_particle & reference_particle);
public:
    FF_drift();
    template <typename T>
    inline static void drift_unit(T & x, T const& xp,
                                  T & y, T const& yp,
                                  T & z, T const& dpop,
                                  double length, double reference_momentum,
                                  double m, double reference_cdt);
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_drift();
};

typedef boost::shared_ptr<FF_drift > FF_drift_sptr;

#include "synergia/utils/invsqrt.h"

template <typename T>
void FF_drift::drift_unit(T & x, T const& xp,
                          T & y, T const& yp,
                          T & cdt, T const& dpop,
                          double length, double reference_momentum,
                          double m, double reference_cdt) {
    T dp = dpop + 1.0;
    T inv_npz = invsqrt(dp * dp - xp * xp - yp * yp);
    T lxpr = xp * length * inv_npz;
    T lypr = yp * length * inv_npz;
    T D2 = lxpr * lxpr + length * length + lypr * lypr;
    T p = dp * reference_momentum;
    T E2 = p * p + m * m;
    //T beta2 = p*p / E2;
    T ibeta2 = E2 / (p * p);
    x += lxpr;
    y += lypr;
    //cdt += sqrt(D2 / beta2) - reference_cdt;
    cdt += ((0.0 < length) - (length < 0.0)) * sqrt(D2 * ibeta2) - reference_cdt;
}
#endif // FF_DRIFT_H
