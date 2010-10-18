#include "reference_particle.h"
#include "synergia/utils/floating_point.h"

Reference_particle::Reference_particle()
{

}

Reference_particle::Reference_particle(int charge, double mass,
        double total_energy) :
    charge(charge), four_momentum(mass, total_energy),
            state(boost::extents[6]), repetition(0), repetition_length(0),
            s(0)
{
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(int charge,
        Four_momentum const & four_momentum_in) :
    charge(charge), four_momentum(four_momentum_in), state(boost::extents[6]),
            repetition(0), repetition_length(0), s(0)
{
    for (int i = 0; i < 6; ++i) {
        state[i] = 0;
    }
}

Reference_particle::Reference_particle(int charge,
        Four_momentum const & four_momentum_in, Const_MArray1d_ref state) :
    charge(charge), four_momentum(four_momentum_in), state(state), repetition(
            0), repetition_length(0), s(0)
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

void
Reference_particle::increment_trajectory(double length)
{
    s += length;
}

void
Reference_particle::start_repetition()
{
    if (repetition_length == 0.0) {
        repetition_length = s;
    }
    if (s > 0.0) {
        repetition += 1;
    }
    s = 0.0;
}

void
Reference_particle::set_trajectory(int repetition, double repetition_length,
        double s)
{
    this->repetition = repetition;
    this->repetition_length = repetition_length;
    this->s = s;
}

int
Reference_particle::get_charge() const
{
    return charge;
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

double
Reference_particle::get_trajectory_length() const
{
    return repetition * repetition_length + s;
}

double
Reference_particle::get_s() const
{
    return s;
}

int
Reference_particle::get_repetition() const
{
    return repetition;
}

double
Reference_particle::get_repetition_length() const
{
    return repetition_length;
}

bool
Reference_particle::equal(Reference_particle const& reference_particle,
        double tolerance) const
{
    if (charge != reference_particle.get_charge()) {
        return false;
    }
    if (!four_momentum.equal(reference_particle.get_four_momentum(), tolerance)) {
        return false;
    }
    for (int i = 0; i < 6; ++i) {
        if (!floating_point_equal(state[i], reference_particle.get_state()[i],
                tolerance)) {
            return false;
        }
    }
    return true;
}
