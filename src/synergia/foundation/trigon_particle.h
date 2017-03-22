#ifndef TRIGON_PARTICLE_H
#define TRIGON_PARTICLE_H

#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/trigon.h"
#include "synergia/utils/multi_array_typedefs.h"

#include <array>
#include <complex>

template <unsigned int Order>
class Trigon_particle
{
public:
    typedef Trigon<double, Order, 6> Component_t;
    typedef Trigon<std::complex<double>, Order, 6> Complex_component_t;
    typedef std::array<Component_t, 6> State_t;

    Trigon_particle(Reference_particle const& reference_particle)
        : reference_particle(reference_particle)
        , state({ { Component_t(reference_particle.get_state()[0], 0),
                    Component_t(reference_particle.get_state()[1], 1),
                    Component_t(reference_particle.get_state()[2], 2),
                    Component_t(reference_particle.get_state()[3], 3),
                    Component_t(reference_particle.get_state()[4], 4),
                    Component_t(reference_particle.get_state()[5], 5) } })
    {
    }

    Reference_particle& get_reference_particle() { return reference_particle; }

    Reference_particle const& get_reference_particle() const
    {
        return reference_particle;
    }

    State_t& get_state() { return state; }

    State_t const& get_state() const { return state; }

    MArray2d get_jacobian() const
    {
        MArray2d retval(boost::extents[6][6]);
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                retval[i][j] = state[i].template get_subpower<1>().terms[j];
            }
        }
        return retval;
    }

private:
    Reference_particle reference_particle;
    State_t state;
};

typedef Trigon_particle<1> Trigon_particle_t;
#endif // TRIGON_PARTICLE_H
