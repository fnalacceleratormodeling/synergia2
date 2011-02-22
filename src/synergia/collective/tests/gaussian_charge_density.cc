#include "gaussian_charge_density.h"

#include <gsl/gsl_sf_erf.h>
#include <cmath>
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;

double
gaussian_charge_density(double Q, double r2, double sigma)
{
    return Q / std::pow(sigma * std::sqrt(2 * pi), 3) * exp(-r2
            / (2 * sigma * sigma));
}

double
gaussian_electric_potential(double Q, double r, double sigma)
{
    return Q / (4.0 * pi * epsilon0 * r) * gsl_sf_erf(r / (std::sqrt(2.0)
            * sigma));
}

double
gaussian_electric_field_component(double Q, double r, double sigma, double x)
{
    return Q * x / (4.0 * pi * epsilon0 * r * r * r) * (std::sqrt(2.0) * r
            / (std::sqrt(pi) * sigma) * exp(-r * r / (2 * sigma * sigma))
            - gsl_sf_erf(r / (std::sqrt(2.0) * sigma)));
}
