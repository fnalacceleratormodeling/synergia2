#ifndef GAUSSIAN_CHARGE_DENSITY_H_
#define GAUSSIAN_CHARGE_DENSITY_H_

// Spherical Gaussian charge density [C/m^3]
double
gaussian_charge_density(double Q, double r2, double sigma);

double
gaussian_charge_density2(double Q, double x, double y, double z, double sigma,
        double sigmaz);

double
gaussian_charge_density3(double Q, double x, double y, double z, double sigmax,
        double sigmay, double sigmaz);

// Electric potential due to spherical Gaussian charge density [V]
double
gaussian_electric_potential(double Q, double r, double sigma);

// Electric field component due to spherical Gaussian charge density [V]
// Since all components are the same, we don't pretend to distinguish
// between them.
double
gaussian_electric_field_component(double Q, double r, double sigma, double x);

double
gaussian_electric_force_component(double q, double Q, double r, double x,
        double sigma);

double
gaussian_electric_field_component2(double Q, double x, double y, double z,
        double sigmax, double sigmay, double sigmaz, double var);

double
gaussian_electric_force_component2(double q, double Q, double x, double y,
        double z, double sigma, double sigmaz, double var);

#endif /* GAUSSIAN_CHARGE_DENSITY_H_ */
