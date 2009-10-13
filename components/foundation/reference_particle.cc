#include "reference_particle.h"

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

Four_momentum const &
Reference_particle::get_four_momentum()
{
    return four_momentum;
}

Const_MArray1d_ref
Reference_particle::get_state()
{
    return state;
}
