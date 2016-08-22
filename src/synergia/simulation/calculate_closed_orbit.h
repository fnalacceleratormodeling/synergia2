#ifndef CALCULATE_CLOSED_ORBIT_H
#define CALCULATE_CLOSED_ORBIT_H

#include "synergia/lattice/lattice.h"

#define DEFAULT_CLOSED_ORBIT_TOLERANCE 1.0e-13

MArray1d
calculate_closed_orbit(Lattice_sptr const lattice_sptr, const double dpp=0.0, const double tolerance=DEFAULT_CLOSED_ORBIT_TOLERANCE);

#endif // CALCULATE_CLOSED_ORBIT_H
