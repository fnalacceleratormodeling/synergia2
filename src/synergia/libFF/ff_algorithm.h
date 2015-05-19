#ifndef FF_RTHICK_H
#define FF_RTHICK_H

#include "ff_drift.h"

class FF_algorithm
{
public:
    template <typename T, void(kf)(T const & x, T & xp, T const & y, T & yp, double const * kL)>
    inline static void yoshida(T & x, T & xp,
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

        const double k1 = d1 * step_strength;
        const double k2 = d2 * step_strength;
        const double k3 = d3 * step_strength;

        for(int i = 0; i < steps; ++i) 
        {
            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c1 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d1 * step_strength );
            kf( x, xp, y, yp, &k1 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c2 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d2 * step_strength );
            kf( x, xp, y, yp, &k2 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c3 * step_length, reference_momentum,
                       m, substep_reference_cdt);

            //kf( x, xp, y, yp, d3 * step_strength );
            kf( x, xp, y, yp, &k3 );

            FF_drift::drift_unit(x, xp, y, yp, cdt, dpop, c4 * step_length, reference_momentum,
                       m, substep_reference_cdt);
        }
    }

};

#endif // FF_QUADRUPOLE_H
