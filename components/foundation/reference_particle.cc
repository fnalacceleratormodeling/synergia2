#include "reference_particle.h"

Reference_particle::Reference_particle(double total_energy, double units[6]) {
	this->total_energy = total_energy;
	for(int i=0; i<6; ++i){
		this->units[i] = units[i];
		this->state[i] = 0;
	}
}

Reference_particle::Reference_particle(double total_energy, double units[6],
		double state[6]) {
	this->total_energy = total_energy;
	for(int i=0; i<6; ++i){
		this->units[i] = units[i];
		this->state[i] = state[i];
	}
}

void
Reference_particle::set_total_energy(double total_energy)
{
	this->total_energy = total_energy;
}

void
Reference_particle::set_state(double state[6])
{
	for(int i=0; i<6; ++i){
		this->state[i] = state[i];
	}
}

double
Reference_particle::get_total_energy()
{
	return total_energy;
}

double *
Reference_particle::get_units()
{
	return units;
}

double *
Reference_particle::get_state()
{
	return state;
}
