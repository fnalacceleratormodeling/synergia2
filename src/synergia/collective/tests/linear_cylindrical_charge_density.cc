#include "linear_cylindrical_charge_density.h"
#include "synergia/foundation/physical_constants.h"
#include <cmath>

double
linear_cylindrical_charge_density(double lambda, double r, double r0)
{
    double retval;
    if (r < r0) {
        retval = 3 * lambda / (mconstants::pi * r0 * r0) * (1.0 - r / r0);
    } else {
        retval = 0.0;
    }
    return retval;
}

// Electric potential due to cylindrical linear charge density [V]
double
linear_cylindrical_electric_potential(double lambda, double r, double r0)
{
    double retval;
    if (r < r0) {
        retval = lambda / (mconstants::pi * pconstants::epsilon0) * (-1.0 / (12
                * r0 * r0 * r0) * (9 * r * r * r0 - 4 * r * r * r) + 5 / 12.0
                - 0.5 * std::log(r0));
    } else {
        retval = -lambda / (2 * mconstants::pi * pconstants::epsilon0)
                * std::log(r);
    }
    return retval;
}

// Transverse electric field component due to cylindrical linear
// charge density [V]
// Since both x- and y-components are the same, we don't pretend to distinguish
// between them.
double
linear_cylindrical_electric_field_component(double lambda, double r, double r0,
        double x)
{
    return 0.0;
}

