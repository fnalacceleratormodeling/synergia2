#include "four_momentum.h"
#include <cmath>

void
Four_momentum::update_from_gamma()
{
	energy= gamma*mass;
	beta = sqrt(1.0-1.0/(gamma*gamma));
	momentum = gamma*beta*mass;
}

Four_momentum::Four_momentum(double mass)
{
    this->mass = mass;
    gamma = 1.0;
    update_from_gamma();
}

Four_momentum::Four_momentum(double mass, double total_energy)
{
    this->mass = mass;
    set_total_energy(total_energy);
}

void
Four_momentum::set_total_energy(double total_energy)
{
	gamma = total_energy/mass;
	update_from_gamma();
}

void
Four_momentum::set_kinetic_energy(double kinetic_energy)
{
	gamma = (mass + kinetic_energy)/mass;
	update_from_gamma();
}

void
Four_momentum::set_momentum(double momentum)
{
	double r2 = momentum*momentum/(mass*mass);
	beta = sqrt(r2/(1+r2));
	set_beta(beta);
}

void
Four_momentum::set_gamma(double gamma)
{
	this->gamma = gamma;
	update_from_gamma();
}

void
Four_momentum::set_beta(double beta)
{
	gamma = 1.0/sqrt(1.0-beta*beta);
	update_from_gamma();
}

double
Four_momentum::get_mass() const
{
	return mass;
}

double
Four_momentum::get_total_energy() const
{
	return energy;
}

double Four_momentum::get_kinetic_energy() const
{
	return energy-mass;
}

double
Four_momentum::get_momentum() const
{
	return momentum;
}

double
Four_momentum::get_gamma() const
{
	return gamma;
}

double
Four_momentum::get_beta() const
{
	return beta;
}
