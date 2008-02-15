#include "constraints.h"
#include <cmath>

const double pi = 3.141592653589793238462643;

void
apply_longitudinal_periodicity(Macro_bunch_store &mbs)
{
    Array_1d<double> z = mbs.local_particles.slice(vector2(Range(4),Range()));
    for (Array_1d<double>::Iterator it = z.begin();
            it != z.end();
            ++it) {
        double tmp = *it + pi;
        if (tmp > 0) {
            *it = fmod(tmp, 2.0 * pi) - pi;
        } else {
            *it = fmod(tmp, 2.0 * pi) + pi;
        }
    }
}

void 
apply_circular_aperture(Macro_bunch_store &mbs)
{
}

