#include "four_momentum.h"
#include <cmath>

void
Four_momentum::update_from_gamma()
{
	energy= gamma*mass;
	beta = 1.0/sqrt(1.0-1.0/(gamma*gamma));
	momentum = gamma*beta*mass;
}

Four_momentum::Four_momentum(double mass, double total_energy)
{
	this->mass = mass;
	set_total_energy(0.0);
}

void
Four_momentum::set_total_energy(double total_energy)
{
	gamma = energy/mass;
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
Four_momentum::get_mass()
{
	return mass;
}

double
Four_momentum::get_total_energy()
{
	return energy;
}

double Four_momentum::get_kinetic_energy()
{
	return energy-mass;
}

double
Four_momentum::get_momentum()
{
	return momentum;
}

double
Four_momentum::get_gamma()
{
	return gamma;
}

double
Four_momentum::get_beta()
{
	return beta;
}
