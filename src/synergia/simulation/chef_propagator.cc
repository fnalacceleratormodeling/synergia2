#include "chef_propagator.h"
#include "synergia/lattice/chef_utils.h"

Chef_propagator::Chef_propagator(Chef_elements const& chef_elements) :
    chef_elements(chef_elements)
{
}

// jfa: This routine is incorrect when passing through an accelerating element.
// Please fix it.
void
Chef_propagator::apply(Bunch & bunch)
{
    std::vector<double > u = chef_unit_conversion(
            bunch.get_reference_particle());
    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    Vector chef_state(6);
    Particle particle(reference_particle_to_chef_particle(
            bunch.get_reference_particle()));
    for (int part = 0; part < local_num; ++part) {
        for (int synergia_index = 0; synergia_index < 6; ++synergia_index) {
            int chef_idx = get_chef_index(synergia_index);
            chef_state[chef_idx] = particles[part][synergia_index]
                    / u[synergia_index];
        }
        particle.State() = chef_state;
        for (Chef_elements::iterator it = chef_elements.begin(); it
                != chef_elements.end(); ++it) {
            (*it)->propagate(particle);
        }
        chef_state = particle.State();
        for (int synergia_index = 0; synergia_index < 6; ++synergia_index) {
            int chef_idx = get_chef_index(synergia_index);
            particles[part][synergia_index] = chef_state[chef_idx]
                    * u[synergia_index];
        }
    }
}
