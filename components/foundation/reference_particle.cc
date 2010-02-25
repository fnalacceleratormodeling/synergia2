#include "reference_particle.h"

Reference_particle::Reference_particle(double mass, double total_energy) :
    four_momentum(mass, total_energy), state(boost::extents[6])
{
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(Four_momentum const & four_momentum_in) :
    four_momentum(four_momentum_in), state(boost::extents[6])
{
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(Four_momentum const & four_momentum_in,
        Const_MArray1d_ref state) :
    four_momentum(four_momentum_in), state(state)
{
}

void
Reference_particle::set_four_momentum(Four_momentum const & four_momentum)
{
    this->four_momentum = four_momentum;
}

void
Reference_particle::set_state(Const_MArray1d_ref state)
{
    this->state = state;
}

void
Reference_particle::set_total_energy(double total_energy)
{
    four_momentum.set_total_energy(total_energy);
}

Four_momentum const &
Reference_particle::get_four_momentum() const
{
    return four_momentum;
}

Const_MArray1d_ref
Reference_particle::get_state() const
{
    return state;
}

double
Reference_particle::get_gamma() const
{
    return four_momentum.get_gamma();
}

double
Reference_particle::get_beta() const
{
    return four_momentum.get_beta();
}

double
Reference_particle::get_momentum() const
{
    return four_momentum.get_momentum();
}

double
Reference_particle::get_total_energy() const
{
    return four_momentum.get_total_energy();
}

bool
Reference_particle::equal(Reference_particle const& reference_particle,
        double tolerance) const
{
    if (!four_momentum.equal(reference_particle.get_four_momentum(), tolerance)) {
        return false;
    }
    for (int i = 0; i < 6; ++i) {
        if (std::abs(state[i]) < tolerance) {
            if (std::abs(state[i] - reference_particle.get_state()[i])
                    > tolerance) {
                return false;
            }
        } else {
            if (std::abs((state[i] - reference_particle.get_state()[i])
                    / state[i]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}
