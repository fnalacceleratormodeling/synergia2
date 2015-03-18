#ifndef FF_SBEND_H
#define FF_SBEND_H

#include "ff_element.h"

class FF_sbend : public FF_element
{
public:
    FF_sbend();

    template <typename T>
    inline static void sbend_unit(T & x, T & xp,
                                  T & y, T & yp,
                                  T & cdt, T const& dpop,
                                  double length, 
                                  double cos_angle, 
                                  double sin_angle,
                                  double reference_momentum,
                                  double m, 
                                  double reference_brho );

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_sbend();
};

typedef boost::shared_ptr<FF_sbend > FF_sbend_sptr;

template <typename T>
inline void FF_sbend::sbend_unit(T & x, T & xp,
                                 T & y, T & yp,
                                 T & cdt, T const& dpop,
                                 double length, 
                                 double cos_angle, 
                                 double sin_angle,
                                 double reference_momentum,
                                 double m, 
                                 double reference_brho )
{
    T rho   = reference_brho * (dpop + 1.0);
    T p     = reference_momentum * (dpop + 1.0);
    T E     = sqrt(p * p + m * m);
    T igamma= m / E;
    T ibeta = E / p;

    T m00 = cos_angle;
    T m03 = sin_angle * rho;
    T m05 = rho * (1.0 - cos_angle) * ibeta;

    T m14 = length;

    T m20 = (-sin_angle) * ibeta;
    T m23 = -m05;
    T m25 = length * ibeta * (igamma * igamma) - (length - sin_angle * rho) * ibeta;

    T m30 = (-sin_angle) / rho;
    T m33 = m00;
    T m35 = -m20;

    T x_ = x;

    x = m00 * x_ + m03 * xp + m05 * dpop;
    y = m14 * yp;
    cdt = m20 * x_ + m23 * xp + m25 * dpop;
    xp  = m30 * x_ + m33 * xp + m35 * dpop;
}

#endif // FF_SBEND_H
