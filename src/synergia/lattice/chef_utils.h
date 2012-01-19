#ifndef CHEF_UTILS_H_
#define CHEF_UTILS_H_

#include <beamline/beamline.h>
#include <beamline/Particle.h>
#include <beamline/JetParticle.h>
#include "synergia/foundation/reference_particle.h"

void
print_chef_beamline(BmlPtr beamline_sptr);

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

Reference_particle
chef_particle_to_reference_particle(Particle const& chef_particle);

JetParticle
reference_particle_to_chef_jet_particle(
        Reference_particle const& reference_particle, int map_order);

Reference_particle
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
inline
int
get_synergia_index(int chef_index)
{
  return 2*(chef_index%3) + chef_index%3;
}
#endif /* CHEF_UTILS_H_ */
