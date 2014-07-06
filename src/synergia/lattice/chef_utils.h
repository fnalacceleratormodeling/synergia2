#ifndef CHEF_UTILS_H_
#define CHEF_UTILS_H_

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <beamline/beamline.h>
#include <beamline/Particle.h>
#include <beamline/JetParticle.h>
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

#include "synergia/foundation/reference_particle.h"

std::string
chef_element_as_string(ElmPtr element_sptr);

std::string
chef_beamline_as_string(BmlPtr beamline_sptr);

void
print_chef_element(ElmPtr element_sptr);

void
print_chef_beamline(BmlPtr beamline_sptr);

std::string
full_chef_beamline_as_string(BmlPtr beamline_sptr);

void
print_full_chef_beamline(BmlPtr beamline_sptr);

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

Reference_particle
chef_particle_to_reference_particle(Particle const& chef_particle);

void
ensure_jet_environment(int map_order);

JetParticle
reference_particle_to_chef_jet_particle(
        Reference_particle const& reference_particle, int map_order);

Reference_particle
propagate_reference_particle(Reference_particle const& reference_particle,
        BmlPtr beamline_sptr);

Particle
get_closed_orbit_particle(Particle util_part, BmlPtr beamline_sptr, double dpop);

/// units conversion
/// X_synergia = U X_chef
/// where U = diag(u[0],u[1],u[2],u[3],u[4],u[5])
std::vector<double >
chef_unit_conversion(Reference_particle const& reference_particle);

inline
int
get_chef_index(int synergia_index)
{
    const int chef_index[] = {0, 3, 1, 4, 2, 5};
    return chef_index[synergia_index];
}

inline
int
get_synergia_index(int chef_index)
{
  const int synergia_index[] = {0, 2, 4, 1, 3, 5};
  return synergia_index[chef_index];
}
#endif /* CHEF_UTILS_H_ */
