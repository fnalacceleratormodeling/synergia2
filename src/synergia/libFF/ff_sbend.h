#ifndef FF_SBEND_H
#define FF_SBEND_H

#include "ff_element.h"
#include "ff_algorithm.h"
#include "synergia/utils/invsqrt.h"

#include <iostream>

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

    template <typename T>
    inline static void sbend_unit2(T & x, T & xp,
                                   T & y, T & yp,
                                   T & cdt, T const& dpop,
                                   double length, 
                                   double angle, 
                                   double strength, 
                                   double reference_momentum,
                                   double m, 
                                   double reference_cdt,
                                   std::complex<double> phase,
                                   std::complex<double> term
                                   );

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    double get_reference_cdt(double length, double strength, double angle, 
                             bool ledge, bool redge,
                             double e1, double e2, double dphi,
                             std::complex<double> const & phase,
                             std::complex<double> const & term,
                             Reference_particle &reference_particle);

    double get_reference_cdt(double length, double angle, double strength,
                             Reference_particle &reference_particle);

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

template <typename T>
inline void FF_sbend::sbend_unit2(T & x, T & xp,
                                 T & y, T & yp,
                                 T & cdt, T const& dpop,
                                 double length, 
                                 double angle, 
                                 double strength, 
                                 double reference_momentum,
                                 double m, 
                                 double reference_cdt,
                                 std::complex<double> phase,
                                 std::complex<double> term
                                 )
{
    typedef std::complex<T> CT;

    T p0 = reference_momentum;
    T p  = reference_momentum * (dpop + 1.0);
    T E0 = sqrt(p0 * p0 + m * m);
    T E  = sqrt(p * p + m * m);

    T igamma = m / E0;
    T ibeta  = invsqrt(1.0 - igamma * igamma);

    T csq = PH_MKS_c * PH_MKS_c * 1e-9;
    T psq = (dpop + 1.0) * (dpop + 1.0);

    T Ef = invsqrt(psq + igamma * igamma * ibeta * ibeta);

    T beta1 = Ef * xp;
    T beta2 = Ef * yp;
    T beta3 = Ef * sqrt( (dpop + 1.0) * (dpop + 1.0) - xp * xp - yp * yp );

    CT ui  = CT(0.0, x);
    CT vui = CT(PH_MKS_c * beta3, PH_MKS_c * beta1);

    T iomega = E / (csq * strength);

    CT bi = CT(0.0, 1.0) * vui * iomega - ui;
    CT bf = bi * phase + term;

    T rho = PH_MKS_c * sqrt( beta1 * beta1 + beta3 * beta3 ) * iomega;

    T dthmphi = asin(bi.real() / rho) - asin(bf.real() / rho);

    CT expf = std::exp( CT(0.0, dthmphi) );
    CT vuf  = vui * expf;
    CT uf   = (ui + bi) * expf - bf;

    T dtheta = dthmphi - angle;
    T ncdt = - PH_MKS_c * dtheta * iomega;

    x    = uf.imag();
    y   += beta2 * ncdt;
    cdt += ncdt - reference_cdt;
    xp   = vuf.imag() / (Ef * PH_MKS_c);
}

#endif // FF_SBEND_H
