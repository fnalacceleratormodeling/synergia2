#include "poly_particle.h"

Poly_particle::Poly_particle(Reference_particle const& reference_particle):
    reference_particle(reference_particle),
    state(6)
{
    for(unsigned int i = 0; i < 6; ++i) {
        state[i].d[i] = 1.0;
    }
}

Reference_particle &
Poly_particle::get_reference_particle()
{
    return reference_particle;
}

Reference_particle const&
Poly_particle::get_reference_particle() const
{
    return reference_particle;
}

Poly_particle::State_t &
Poly_particle::get_state()
{
    return state;
}

Poly_particle::State_t const&
Poly_particle::get_state() const
{
    return state;
}

MArray2d
Poly_particle::get_jacobian() const
{
    MArray2d retval(boost::extents[6][6]);
    for(unsigned int i=0; i<6; ++i) {
        for(unsigned int j=0; j<=i; ++j) {
            retval[i][j] = retval[j][i] = state[i].d[j];
        }
    }
    return retval;
}
