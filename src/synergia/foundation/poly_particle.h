#ifndef POLY_PARTICLE_H
#define POLY_PARTICLE_H

#include "synergia/foundation/poly.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/utils/multi_array_typedefs.h"

#include <vector>

class Poly_particle
{
public:
    typedef Poly Component_t;
    typedef std::vector<Component_t> State_t;

    Poly_particle(Reference_particle const& reference_particle);
    Reference_particle & get_reference_particle();
    Reference_particle const& get_reference_particle() const;
    State_t & get_state();
    State_t const& get_state() const;
    MArray2d get_jacobian() const;

private:
    Reference_particle reference_particle;
    State_t state;
};

#endif // POLY_PARTICLE_H
