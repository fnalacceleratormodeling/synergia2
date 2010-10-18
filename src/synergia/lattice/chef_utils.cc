#include "chef_utils.h"
#include <iostream>
#include "synergia/utils/floating_point.h"

void
print_chef_beamline(BmlPtr beamline_sptr)
{
    for (beamline::const_iterator it = beamline_sptr->begin(); it
            != beamline_sptr->end(); ++it) {
        std::cout << (*it)->Name() << "(" << (*it)->Type() << "): Length="
                << (*it)->Length() << ", Strength=" << (*it)->Strength()
                << std::endl;
    }
}

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle)
{
    // n.b. We don't deal with negatively charged particles!!!!
    // n.b. We also don't deal with non-zero state.
    double mass = reference_particle.get_four_momentum().get_mass();
    double momentum = reference_particle.get_momentum();
    const double mass_tolerance = 1.0e-5;
    if (floating_point_equal(mass, PH_NORM_mp, mass_tolerance)) {
        Proton proton;
        proton.SetReferenceMomentum(momentum);
        proton.setStateToZero();
        return proton;
    } else {
        if (floating_point_equal(mass, PH_NORM_me, mass_tolerance)) {
            Positron positron;
            positron.SetReferenceMomentum(momentum);
            positron.setStateToZero();
            return positron;
        } else {
            if (floating_point_equal(mass, PH_NORM_mmu, mass_tolerance)) {
                AntiMuon antimuon;
                antimuon.SetReferenceMomentum(momentum);
                antimuon.setStateToZero();
                return antimuon;
            } else {
                throw(runtime_error(
                        "reference_particle_to_chef_particle: particle mass not equal to proton, electron or muon mass."));
            }
        }
    }
}

JetParticle
reference_particle_to_chef_jet_particle(
        Reference_particle const& reference_particle, int map_order)
{
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
    return JetParticle(reference_particle_to_chef_particle(reference_particle));
}

void
propagate_reference_particle(Reference_particle const& reference_particle,
        BmlPtr beamline_sptr)
{
    Particle particle(reference_particle_to_chef_particle(reference_particle));
    std::cout << "initial state: " << particle.State() << std::endl;
    beamline_sptr->propagate(particle);
    std::cout << "final state: " << particle.State() << std::endl;
}

std::vector<double >
chef_unit_conversion(Reference_particle const& reference_particle)
{
    std::vector<double > u(6);
    u.at(0) = 1.0;
    u.at(1) = reference_particle.get_momentum();
    u.at(2) = 1.0;
    u.at(3) = reference_particle.get_momentum();
    u.at(4) = 1.0;
    u.at(5) = reference_particle.get_momentum();
    return u;
}

