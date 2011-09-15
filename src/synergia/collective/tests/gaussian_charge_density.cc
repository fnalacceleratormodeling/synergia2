#include "gaussian_charge_density.h"

#include <iostream>
#include <gsl/gsl_sf_erf.h>
#include <cmath>
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;

// sigmax = sigmay = sigmaz
double
gaussian_charge_density(double Q, double r2, double sigma)
{
    return Q / std::pow(sigma * std::sqrt(2 * pi), 3) * exp(-r2 / (2 * sigma
            * sigma));
}

// sigmax = sigmay != sigmaz
double
gaussian_charge_density2(double Q, double x, double y, double z, double sigma,
        double sigmaz)
{
    return Q / std::pow(sigma * std::sqrt(2 * pi), 2) * exp(-(x * x + y * y)
            / (2 * sigma * sigma)) * exp (- z * z / (2 * sigmaz * sigmaz))
            / (std::sqrt(2 * pi) * sigmaz);
}

// sigmax != sigmay != sigmaz
double
gaussian_charge_density3(double Q, double x, double y, double z, double sigmax,
        double sigmay, double sigmaz)
{
    return Q / (std::pow(std::sqrt(2 * pi), 3) * sigmax * sigmay * sigmaz)
            * exp(-(x * x) / (2 * sigmax * sigmax) - (y * y) / (2 * sigmay
            * sigmay) - (z * z) / (2 * sigmaz * sigmaz));
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
    return - Q * x / (4.0 * pi * epsilon0 * r * r * r) * (std::sqrt(2.0) * r
            / (std::sqrt(pi) * sigma) * exp(-r * r / (2 * sigma * sigma))
            - gsl_sf_erf(r / (std::sqrt(2.0) * sigma)));
}

double
gaussian_electric_force_component(double q, double Q, double r, double x,
        double sigma)
{
    return - q * Q * x / (4.0 * pi * epsilon0 * r * r * r) * (std::sqrt(2.0)
            * r / (std::sqrt(pi) * sigma) * exp(-r * r / (2 * sigma * sigma))
            - gsl_sf_erf(r / (std::sqrt(2.0) * sigma)));
}

double
gaussian_electric_field_component2(double Q, double x, double y, double z,
        double sigmax, double sigmay, double sigmaz, double var)
{
    double dt = 0.001;
    int NT = 1000;
    double integral_part = 0.0;
    for (int i = 0; i < NT; ++i) {
        double t = dt * i;
        double lambda = t / (1.0 - t);
        double lambda_square = lambda * lambda;

        double lambda_sigma_x = lambda_square * sigmax * sigmax + 1.0;
        double lambda_sigma_y = lambda_square * sigmay * sigmay + 1.0;
        double lambda_sigma_z = lambda_square * sigmaz * sigmaz + 1.0;

        double exponent_x = - x * x * lambda_square / (2.0 * lambda_sigma_x);
        double exponent_y = - y * y * lambda_square / (2.0 * lambda_sigma_y);
        double exponent_z = - z * z * lambda_square / (2.0 * lambda_sigma_z);

        double lambda_sigma_var;
        if (var == x) {
            lambda_sigma_var = lambda_sigma_x;
        } else if (var == y) {
            lambda_sigma_var = lambda_sigma_y;
        }

        // using the simpson's rule for the numerical integration
        if (i % 2 == 0) {
             integral_part += 2.0 * (lambda_square / lambda_sigma_var) 
                * exp(exponent_x) * exp(exponent_y) * exp(exponent_z) 
                / sqrt(lambda_sigma_x * lambda_sigma_y * lambda_sigma_z) 
                / ((1.0 - t) * (1.0 - t));;
        } else if (i % 2 == 1) {
             integral_part += 4.0 * (lambda_square / lambda_sigma_var) 
                * exp(exponent_x) * exp(exponent_y) * exp(exponent_z) 
                / sqrt(lambda_sigma_x * lambda_sigma_y * lambda_sigma_z) 
                / ((1.0 - t) * (1.0 - t));;
        }
    }
    return Q / (4.0 * pi * epsilon0) * std::sqrt(2.0 / pi) * var *
            integral_part * dt / 3.0;
}

double
gaussian_electric_force_component2(double q, double Q, double x, double y,
        double z, double sigma, double sigmaz, double var)
{
    double dt = 0.0001;
    int NT = 10000;
    double integral_part = 0.0;
    for (int i = 0; i < NT; ++i) {
        double t = dt * i;
        double lambda = t / (1.0 - t);
        double c0 = lambda * lambda * sigmaz * sigmaz + 1.0;
        double c1 = lambda * lambda * sigma * sigma + 1.0;
        double c2 = lambda * lambda / 2.0 / c1;
        double c3 = lambda * lambda / 2.0 / c0;
        double exp1 = exp( - x * x * c2);
        double exp2 = exp( - y * y * c2);
        double exp3 = exp( - z * z * c3);
        if (i == 0) {
            integral_part += c2 * exp1 * exp2 * exp3
            / std::sqrt(c1 * c1 * c0) * dt / (1.0 - t) / (1.0 - t);
        } else {
            integral_part += (2.0 * c2) * exp1 * exp2 * exp3
                    / std::sqrt(c1 * c1 * c0) * dt / (1.0 - t) / (1.0 - t);
        }
    }
    return q * Q / (4.0 * pi * epsilon0) * std::sqrt(2.0 / pi) * var *
            integral_part;
}
