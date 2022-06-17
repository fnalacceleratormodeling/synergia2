#ifndef LINEAR_CYLINDRICAL_CHARGE_DENSITY_H_
#define LINEAR_CYLINDRICAL_CHARGE_DENSITY_H_

// Cylindrical linear charge density [C/m^3]
double
linear_cylindrical_charge_density(double lambda, double r, double r0);

// Electric potential due to cylindrical linear charge density [V]
double
linear_cylindrical_electric_potential(double lambda, double r, double r0);

// Transverse electric field component due to cylindrical linear
// charge density [V]
// Since both x- and y-components are the same, we don't pretend to distinguish
// between them.
double
linear_cylindrical_electric_field_component(double lambda, double r, double r0,
        double x);

double
linear_cylindrical_electric_force_component(double q, double lambda, double r,
        double r0, double x);

#endif /* LINEAR_CYLINDRICAL_CHARGE_DENSITY_H_ */
