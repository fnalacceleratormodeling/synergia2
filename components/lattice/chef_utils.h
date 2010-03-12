#ifndef CHEF_UTILS_H_
#define CHEF_UTILS_H_

#include <beamline/beamline.h>
#include <beamline/Particle.h>
#include "components/foundation/reference_particle.h"

void
print_chef_beamline(BmlPtr beamline_sptr);

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

#endif /* CHEF_UTILS_H_ */
