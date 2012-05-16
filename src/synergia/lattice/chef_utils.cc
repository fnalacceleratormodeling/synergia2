#include "chef_utils.h"
#include <iostream>
#include <sstream>
#include "synergia/utils/floating_point.h"

std::string
chef_beamline_as_string(BmlPtr beamline_sptr)
{
    std::stringstream sstream;
    for (beamline::const_iterator it = beamline_sptr->begin(); it
            != beamline_sptr->end(); ++it) {
        sstream << (*it)->Name() << "(" << (*it)->Type() << "): Length="
                << (*it)->Length() << ", Strength=" << (*it)->Strength()
                << std::endl;
    }
    return sstream.str();
}

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
    double mass = reference_particle.get_four_momentum().get_mass();
    double momentum = reference_particle.get_momentum();
    int charge = reference_particle.get_charge();
    const double mass_tolerance = 1.0e-5;
    Vector chef_state(6);
    chef_state[0] = reference_particle.get_state()[0];
    chef_state[1] = reference_particle.get_state()[2];
    chef_state[2] = reference_particle.get_state()[4];
    chef_state[3] = reference_particle.get_state()[1];
    chef_state[4] = reference_particle.get_state()[3];
    chef_state[5] = reference_particle.get_state()[5];
    if (charge == 1) {
        if (floating_point_equal(mass, PH_NORM_mp, mass_tolerance)) {
            Proton proton;
            proton.SetReferenceMomentum(momentum);
            proton.State() = chef_state;
            return proton;
        } else {
            if (floating_point_equal(mass, PH_NORM_me, mass_tolerance)) {
                Positron positron;
                positron.SetReferenceMomentum(momentum);
                positron.State() = chef_state;
                return positron;
            } else {
                if (floating_point_equal(mass, PH_NORM_mmu, mass_tolerance)) {
                    AntiMuon antimuon;
                    antimuon.SetReferenceMomentum(momentum);
                    antimuon.State() = chef_state;
                    return antimuon;
                } else {
                    throw(runtime_error(
                            "reference_particle_to_chef_particle: particle mass not equal to proton, electron or muon mass."));
                }
            }
        }
    } else {
        if (charge == -1) {
            if (floating_point_equal(mass, PH_NORM_mp, mass_tolerance)) {
                AntiProton antiproton;
                antiproton.SetReferenceMomentum(momentum);
                antiproton.State() = chef_state;
                return antiproton;
            } else {
                if (floating_point_equal(mass, PH_NORM_me, mass_tolerance)) {
                    Electron electron;
                    electron.SetReferenceMomentum(momentum);
                    electron.State() = chef_state;
                    return electron;
                } else {
                    if (floating_point_equal(mass, PH_NORM_mmu, mass_tolerance)) {
                        Muon muon;
                        muon.SetReferenceMomentum(momentum);
                        muon.State() = chef_state;
                        return muon;
                    } else {
                        throw(runtime_error(
                                "reference_particle_to_chef_particle: particle mass not equal to proton, electron or muon mass."));
                    }
                }
            }
        } else {
            throw(runtime_error(
                    "reference_particle_to_chef_particle: particle does not have unit charge."));
        }
    }
}

Reference_particle
chef_particle_to_reference_particle(Particle const& chef_particle)
{
    Four_momentum four_momentum(chef_particle.Mass());
    four_momentum.set_momentum(chef_particle.Momentum());
    Reference_particle
            reference_particle(chef_particle.Charge(), four_momentum);
    MArray1d state(boost::extents[6]);
    state[0] = chef_particle.State()[0];
    state[1] = chef_particle.State()[3];
    state[2] = chef_particle.State()[1];
    state[3] = chef_particle.State()[4];
    state[4] = chef_particle.State()[2];
    state[5] = chef_particle.State()[5];
    reference_particle.set_state(state);
    return reference_particle;
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

Reference_particle
propagate_reference_particle(Reference_particle const& reference_particle,
        BmlPtr beamline_sptr)
{
    Particle particle(reference_particle_to_chef_particle(reference_particle));
    beamline_sptr->propagate(particle);
    return chef_particle_to_reference_particle(particle);
}

std::vector<double >
chef_unit_conversion(Reference_particle const& reference_particle)
{
    std::vector<double > u(6);
    u[0] = 1.0;
    u[1] = 1.0;
    u[2] = 1.0;
    u[3] = 1.0;
    u[4] = 1.0;
    u[5] = 1.0;
    return u;
}

