#ifndef CHEF_UTILS_H_
#define CHEF_UTILS_H_

#include <beamline/beamline.h>
#include <beamline/Particle.h>
#include <beamline/JetParticle.h>
#include "components/foundation/reference_particle.h"

void
print_chef_beamline(BmlPtr beamline_sptr);

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

JetParticle
reference_particle_to_chef_jet_particle(
        Reference_particle const& reference_particle, int map_order);

void
propagate_reference_particle(Reference_particle const& reference_particle,
        BmlPtr beamline_sptr);

/// units conversion
/// X_synergia = U X_chef
/// where U = diag(u[0],u[1],u[2],u[3],u[4],u[5])
std::vector<double >
chef_unit_conversion(Reference_particle const& reference_particle);

inline
int
get_chef_index(int synergia_index)
{
    return synergia_index / 2 + 3 * (synergia_index % 2);
}

#endif /* CHEF_UTILS_H_ */
