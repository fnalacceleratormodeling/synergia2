#include "reference_particle.h"

Reference_particle::Reference_particle(double energy) {
	this->energy = energy;
	for(int i=0; i<6; ++i){
		this->state[i] = 0;
	}
}

Reference_particle::Reference_particle(double energy, double state[6]) {
	this->energy = energy;
	for(int i=0; i<6; ++i){
		this->state[i] = state[i];
	}
}

void
Reference_particle::set_energy(double energy)
{
	this->energy = energy;
}

void
Reference_particle::set_state(double state[6])
{
	for(int i=0; i<6; ++i){
		this->state[i] = state[i];
	}
}

double
Reference_particle::get_energy()
{
	return energy;
}

double *
Reference_particle::get_state()
{
	return state;
}

