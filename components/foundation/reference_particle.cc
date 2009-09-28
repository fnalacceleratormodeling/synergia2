#include "reference_particle.h"

Reference_particle::Reference_particle(double total_energy) :
    state(boost::extents[6])
{
    this->total_energy = total_energy;
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(double total_energy,
        Const_MArray1d_ref state) :
    state(state)
{
    this->total_energy = total_energy;
}

void
Reference_particle::set_total_energy(double total_energy)
{
    this->total_energy = total_energy;
}

void
Reference_particle::set_state(Const_MArray1d_ref state)
{
    this->state = state;
}

double
Reference_particle::get_total_energy()
{
    return total_energy;
}

MArray1d_ref
Reference_particle::get_state()
{
    return state;
}
