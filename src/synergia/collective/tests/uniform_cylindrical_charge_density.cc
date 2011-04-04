#include "uniform_cylindrical_charge_density.h"
#include "synergia/foundation/physical_constants.h"
#include <cmath>

double
uniform_cylindrical_charge_density(double lambda, double r, double r0)
{
    double retval;
    if (r < r0) {
        retval = lambda / (mconstants::pi * r0 * r0);
    } else {
        retval = 0.0;
    }
    return retval;
}

// Electric potential due to cylindrical uniform charge density [V]
double
uniform_cylindrical_electric_potential(double lambda, double r, double r0)
{
    double retval;
    if (r < r0) {
        retval = -lambda / (4 * mconstants::pi * pconstants::epsilon0) * (r * r
                / (r0 * r0) + 2 * std::log(r0) - 1);
    } else {
        retval = -lambda / (2 * mconstants::pi * pconstants::epsilon0)
                * std::log(r);
    }
    return retval;
}

// Transverse electric field component due to cylindrical uniform
// charge density [V]
// Since both x- and y-components are the same, we don't pretend to distinguish
// between them.
double
uniform_cylindrical_electric_field_component(double lambda, double r,
        double r0, double x)
{
    double retval;
    if (r < r0) {
        retval = lambda / (2 * mconstants::pi * pconstants::epsilon0 * r0 * r0)
                * x;
    } else {
        retval = lambda / (2 * mconstants::pi * pconstants::epsilon0) * (1 / r)
                * x / r;
    }
    return retval;
}

