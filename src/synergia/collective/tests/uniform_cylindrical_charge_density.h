#ifndef UNIFORM_CYLINDRICAL_CHARGE_DENSITY_H_
#define UNIFORM_CYLINDRICAL_CHARGE_DENSITY_H_

// Cylindrical uniform charge density [C/m^3]
double
uniform_cylindrical_charge_density(double lambda, double r, double r0);

// Electric potential due to cylindrical uniform charge density [V]
double
uniform_cylindrical_electric_potential(double lambda, double r, double r0);

// Transverse electric field component due to cylindrical uniform
// charge density [V]
// Since both x- and y-components are the same, we don't pretend to distinguish
// between them.
double
uniform_cylindrical_electric_field_component(double lambda, double r,
        double r0, double x);

double
uniform_cylindrical_electric_force_component(double q, double lambda, double r,
        double r0, double x);

#endif /* UNIFORM_CYLINDRICAL_CHARGE_DENSITY_H_ */
