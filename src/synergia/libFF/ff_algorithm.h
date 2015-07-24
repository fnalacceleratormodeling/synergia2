#ifndef FF_ALGORITHM_H
#define FF_ALGORITHM_H

#include <math.h>
#include "synergia/utils/invsqrt.h"

class FF_algorithm
{
public:

    template <typename T>
    inline static void drift_unit
      (T & x, T const& xp, T & y, T const& yp, T & cdt, T const& dpop,
       double length, double reference_momentum, double m, double reference_cdt) 
    {
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


    template <typename T>
    inline static void thin_dipole_unit
      (T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        xp += -kL[0];
        yp += -kL[1];
    }

    template <typename T>
    inline static void thin_quadrupole_unit
      (T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        xp += -kL[0] * x - kL[1] * y;
        yp += kL[0] * y - kL[1] * x;
    }

    template <typename T>
    inline static void thin_sextupole_unit
      (T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        xp += -0.5 * kL[0] * (x * x - y * y) - kL[1] * x * y;
        yp += kL[0] * x * y - 0.5 * kL[1] * (x * x - y * y);
    }

    template <typename T>
    inline static void thin_rbend_unit
      (T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        thin_dipole_unit(x, xp, y, yp, kL);
        thin_quadrupole_unit(x, xp, y, yp, kL + 2);
        thin_sextupole_unit(x, xp, y, yp, kL + 4);
    }

    // general thin magnets of n-th order
    // n = 1, dipole; n = 2, quadrupole; n = 3, sextupole, etc.
    template <typename T>
    inline static void thin_magnet_unit
      (T const& x, T& xp, T const& y, T& yp, double const * kL, int n)
    {
        for(int k = 0; k < n; k += 4)
        {
            xp += -kL[0] * (n-k) * pow(x, n-k-1) * pow(y, k) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 2; k < n; k += 4)
        {
            xp += +kL[0] * (n-k) * pow(x, n-k-1) * pow(y, k) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 1; k < n; k += 4)
        {
            xp += -kL[1] * (n-k) * pow(x, n-k-1) * pow(y, k) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 3; k < n; k += 4)
        {
            xp += +kL[1] * (n-k) * pow(x, n-k-1) * pow(y, k) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 4; k <= n; k += 4)
        {
            yp += -kL[0] * k * pow(x, n-k) * pow(y, k-1) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 2; k <= n; k += 4)
        {
            yp += +kL[0] * k * pow(x, n-k) * pow(y, k-1) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 1; k <= n; k += 4)
        {
            yp += +kL[1] * k * pow(x, n-k) * pow(y, k-1) 
                         / (factorial(n-k) * factorial(k));
        }

        for(int k = 3; k <= n; k += 4)
        {
            yp += +kL[1] * k * pow(x, n-k) * pow(y, k-1) 
                         / (factorial(n-k) * factorial(k));
        }
    }

    inline double factorial(int n)
    {
        if (n == 0) return 1.0;

        double r = 1; 
        for(int i = 1; i <= n; ++i) r *= i;
        return r;
    }



    // utility
    inline int full_drifts_per_step(int order)
    { return pow(3.0, (order-2.0)/2.0) * 2; }

    inline int compact_drifts_per_step(int order)
    { return (full_drifts_per_step(order) - 2) / 2 + 2; }

    // general n-th order yoshida
    template <
        typename T,
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int n,
        int components >
    struct yoshida_element
    {
        static void integral ( T & x, T & xp,
                               T & y, T & yp,
                               T & cdt, T const & dpop,
                               double pref, double m, double step_ref_cdt,
                               double step_length, double * step_strength,
                               int steps, double c )
        {
            double substep_ref_cdt = step_ref_cdt / 3.0;

            yoshida_element<T, kf, n-1, components>::integral( x, xp, y, yp, cdt, dpop, pref, m, substep_ref_cdt,
                                         step_length, step_strength, steps, c * x1(n) );

            yoshida_element<T, kf, n-1, components>::integral( x, xp, y, yp, cdt, dpop, pref, m, substep_ref_cdt,
                                         step_length, step_strength, steps, c * x0(n) );

            yoshida_element<T, kf, n-1, components>::integral( x, xp, y, yp, cdt, dpop, pref, m, substep_ref_cdt,
                                         step_length, step_strength, steps, c * x1(n) );
        }

        static double x1(int nn)
        { return 1.0 / (2.0 - pow(2.0, 1.0/(2*nn+1))); }

        static double x0(int nn)
        { return 1.0 - 2.0 * x1(nn); }
    };

#if 0
    template <
        typename T,
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int components >
    struct yoshida_element < T, kf, 0, components >
    {
        static void integral ( T & x, T & xp,
                               T & y, T & yp,
                               T & cdt, T const & dpop,
                               double pref, double m, double step_ref_cdt,
                               double step_length, double * step_strength,
                               int steps, double c )
        {
            double substep_ref_cdt = step_ref_cdt / 2.0;
            double kl[components * 2];

            for (int i = 0; i < components * 2; ++i)
                kl[i] = step_strength[i] * c;

            drift_unit( x, xp, y, yp, cdt, dpop, 0.5 * c * step_length, pref, m, substep_ref_cdt );

            kf( x, xp, y, yp, kl );

            drift_unit( x, xp, y, yp, cdt, dpop, 0.5 * c * step_length, pref, m, substep_ref_cdt );
        }
    };
#endif

    template <
        typename T,
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int components >
    struct yoshida_element < T, kf, 0, components >
    {
        static void integral ( T & x, T & xp,
                               T & y, T & yp,
                               T & cdt, T const & dpop,
                               double pref, double m, double step_ref_cdt,
                               double step_length, double * step_strength,
                               int steps, double c )
        {
            double substep_ref_cdt = step_ref_cdt;
            double kl[components * 2];

            for (int i = 0; i < components * 2; ++i)
                kl[i] = step_strength[i] * c * 0.5;

            kf( x, xp, y, yp, kl );

            drift_unit( x, xp, y, yp, cdt, dpop, c * step_length, pref, m, substep_ref_cdt );

            kf( x, xp, y, yp, kl );
        }
    };


    template <
        typename T,
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int order,
        int components >
    inline static void yoshida( T & x, T & xp,
                                T & y, T & yp,
                                T & cdt, T const & dpop,
                                double pref, double m, double step_ref_cdt,
                                double step_length, double * step_strength, int steps )
    {
        const int n = (order - 2) / 2;

        for(int i = 0; i < steps; ++i)
        {
            yoshida_element<T, kf, n, components>::integral( x, xp, y, yp, cdt, dpop, pref, m, step_ref_cdt,
                                                   step_length, step_strength, steps, 1.0 );
        }
    }




    // hardwired 2nd order yoshida
    template <
        typename T, 
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int components >
    inline static void yoshida2(T & x, T & xp,
                                T & y, T & yp,
                                T & cdt, T const& dpop,
                                double reference_momentum,
                                double m, double substep_reference_cdt,
                                double step_length, double * step_strength,
                                int steps)
    {
        for(int i = 0; i < steps; ++i) 
        {
            drift_unit(x, xp, y, yp, cdt, dpop, 0.5 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, step_strength );

            drift_unit(x, xp, y, yp, cdt, dpop, 0.5 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }


    // hardwired 4th order yoshida
    template <
        typename T, 
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int components >
    inline static void yoshida4(T & x, T & xp,
                                T & y, T & yp,
                                T & cdt, T const& dpop,
                                double reference_momentum,
                                double m, double substep_reference_cdt,
                                double step_length, double * step_strength,
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

        double k1[components * 2], k2[components * 2], k3[components * 2];

        for (int i=0; i<components; ++i)
        {
            k1[i*2+0] = d1 * step_strength[i*2+0];
            k1[i*2+1] = d1 * step_strength[i*2+1];

            k2[i*2+0] = d2 * step_strength[i*2+0];
            k2[i*2+1] = d2 * step_strength[i*2+1];

            k3[i*2+0] = d3 * step_strength[i*2+0];
            k3[i*2+1] = d3 * step_strength[i*2+1];
        }

        for(int i = 0; i < steps; ++i) 
        {
            drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d1 * step_strength );
            kf( x, xp, y, yp, k1 );

            drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d2 * step_strength );
            kf( x, xp, y, yp, k2 );

            drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d3 * step_strength );
            kf( x, xp, y, yp, k3 );

            drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }


    // hardwired 6th order yoshida
    template <
        typename T, 
        void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL),
        int components >
    inline static void yoshida6(T & x, T & xp,
                                T & y, T & yp,
                                T & cdt, T const& dpop,
                                double reference_momentum,
                                double m, double substep_reference_cdt,
                                double step_length, double * step_strength,
                                int steps)
    {
        // see yoshida4.py for formulas
        const double c1 = 0.79361246386112147294625603763;
        const double c2 = -0.206276584816439780698721319354;
        const double c3 = -0.118008867881292655922412133144;
        const double c4 = 0.236949573653050744373598734219;

        const double d1 = 1.58722492772224294589251207526;
        const double d2 = -1.99977809735512250728995471397;
        const double d3 = -1.82324266348482825773733634154;
        const double d4 = 2.29714181079092974648453380998;

        double k1[components * 2], k2[components * 2], k3[components * 2], k4[components * 2];

        for (int i=0; i<components; ++i)
        {
            k1[i*2+0] = d1 * step_strength[i*2+0];
            k1[i*2+1] = d1 * step_strength[i*2+1];

            k2[i*2+0] = d2 * step_strength[i*2+0];
            k2[i*2+1] = d2 * step_strength[i*2+1];

            k3[i*2+0] = d3 * step_strength[i*2+0];
            k3[i*2+1] = d3 * step_strength[i*2+1];

            k4[i*2+0] = d4 * step_strength[i*2+0];
            k4[i*2+1] = d4 * step_strength[i*2+1];
        }

        for(int i = 0; i < steps; ++i) 
        {
            drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k2 );

            drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k3 );

            drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k4 );

            drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k3 );

            drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k2 );

            drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }

};

#endif // FF_ALGORITHM_H
