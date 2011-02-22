#ifndef GAUSSIAN_CHARGE_DENSITY_H_
#define GAUSSIAN_CHARGE_DENSITY_H_

// Spherical Gaussian charge density [C/m^3]
double
gaussian_charge_density(double Q, double r2, double sigma);

// Electric potential due to spherical Gaussian charge density [V]
double
gaussian_electric_potential(double Q, double r, double sigma);

// Electric field component due to spherical Gaussian charge density [V]
// Since all components are the same, we don't pretend to distinguish
// between them.
double
gaussian_electric_field_component(double Q, double r, double sigma, double x);

#endif /* GAUSSIAN_CHARGE_DENSITY_H_ */
