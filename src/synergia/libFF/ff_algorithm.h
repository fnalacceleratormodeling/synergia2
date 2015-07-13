#ifndef FF_RTHICK_H
#define FF_RTHICK_H

#include "ff_drift.h"

class FF_algorithm
{
public:

    template <typename T>
    inline static void thin_dipole_unit(T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        xp += -kL[0];
        yp += -kL[1];
    }

    template <typename T>
    inline static void thin_quadrupole_unit(T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        xp += -kL[0] * x - kL[1] * y;
        yp += kL[0] * y - kL[1] * x;
    }

    template <typename T>
    inline static void thin_sextupole_unit(T const& x, T& xp, T const& y, T& yp, double const * kL) 
    {
        xp += -0.5 * kL[0] * (x * x - y * y) - kL[1] * x * y;
        yp += kL[0] * x * y - 0.5 * kL[1] * (x * x - y * y);
    }

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
            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, 0.5 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, step_strength );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, 0.5 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }


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
            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d1 * step_strength );
            kf( x, xp, y, yp, k1 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d2 * step_strength );
            kf( x, xp, y, yp, k2 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d3 * step_strength );
            kf( x, xp, y, yp, k3 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }


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
            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k2 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k3 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k4 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k3 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k2 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            kf( x, xp, y, yp, k1 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }

};

#endif // FF_ALGORITHM_H
