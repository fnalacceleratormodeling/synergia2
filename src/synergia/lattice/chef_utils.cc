#include "chef_utils.h"
#include <iostream>
#include <sstream>
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/math_constants.h"
#include <beamline/beamline.h>
#include <beamline/bmlnElmnt.h>
#include <beamline/CF_sbend.h>
#include <beamline/CF_rbend.h>
#include <beamline/rfcavity.h>
#include <physics_toolkit/ClosedOrbitSage.h>

std::string
chef_beamline_as_string(BmlPtr beamline_sptr)
{
    std::stringstream sstream;
    for (beamline::const_iterator it = beamline_sptr->begin();
            it != beamline_sptr->end(); ++it) {
        sstream << chef_element_as_string(*it);
    }
    return sstream.str();
}

std::string
chef_element_as_string(ElmPtr element_sptr)
{
    std::stringstream sstream;
    sstream << (element_sptr)->Name() << "(" << (element_sptr)->Type() << "): Length="
            << (element_sptr)->Length() << ", Strength=" << (element_sptr)->Strength();
    if ( (std::strcmp((element_sptr)->Type(),"CF_rbend") == 0) ) {
        sstream << ", Quadrupole="
            << boost::static_pointer_cast<CF_rbend>(element_sptr)->getQuadrupole()
                << ", Sextupole="
                << boost::static_pointer_cast<CF_rbend>(element_sptr)->getSextupole();
    } else if ( (std::strcmp((element_sptr)->Type(), "CF_sbend") == 0)) {
        sstream << ", Quadrupole="
            << boost::static_pointer_cast<CF_sbend>(element_sptr)->getQuadrupole()
                << ", Sextupole="
                << boost::static_pointer_cast<CF_sbend>(element_sptr)->getSextupole();
    } else if ( std::strcmp((element_sptr)->Type(), "rfcavity") == 0 ) {
        sstream << ", Freq="
                << boost::static_pointer_cast<rfcavity>(element_sptr)->getRadialFrequency()/(2.0*mconstants::pi) <<
                ", Harmon="
                << boost::static_pointer_cast<rfcavity>(element_sptr)->getHarmonicNumber();
    } else if ( std::strcmp((element_sptr)->Type(), "thinrfcavity") == 0) {
        sstream << ", Freq="
                << boost::static_pointer_cast<thinrfcavity>(element_sptr)->getRadialFrequency()/(2.0*mconstants::pi) <<
                ", Harmon="
                << boost::static_pointer_cast<thinrfcavity>(element_sptr)->getHarmonicNumber();
    }
    sstream << std::endl;
    return sstream.str();
}

void
print_chef_beamline(BmlPtr beamline_sptr)
{
    std::cout << chef_beamline_as_string(beamline_sptr) << std::endl;
    std::cout.flush();
}

void
print_chef_element(ElmPtr element_sptr)
{
    std::cout << chef_element_as_string(element_sptr) << std::endl;
    std::cout.flush();
}

std::string
full_chef_beamline_as_string(BmlPtr beamline_sptr)
{
    std::stringstream sstream;
    for (beamline::deep_iterator it = beamline_sptr->deep_begin();
            it != beamline_sptr->deep_end(); ++it) {
        sstream << (*it)->Name() << "(" << (*it)->Type() << "): Length="
                << (*it)->Length() << ", Strength=" << (*it)->Strength();
        if ( (std::strcmp((*it)->Type(), "CF_rbend") == 0) ) {
            sstream << ", Quadrupole="
                << boost::static_pointer_cast<CF_rbend>(*it)->getQuadrupole()
                    << ", Sextupole="
                    << boost::static_pointer_cast<CF_rbend>(*it)->getSextupole();
        } else if ( (std::strcmp((*it)->Type(), "CF_sbend") == 0)) {
            sstream << ", Quadrupole="
                << boost::static_pointer_cast<CF_sbend>(*it)->getQuadrupole()
                    << ", Sextupole="
                    << boost::static_pointer_cast<CF_sbend>(*it)->getSextupole();
        } else if ( std::strcmp((*it)->Type(), "rfcavity") == 0 ) {
            sstream << ", Freq="
                    << boost::static_pointer_cast<rfcavity>(*it)->getRadialFrequency()/(2.0*mconstants::pi) <<
                    ", Harmon="
                    << boost::static_pointer_cast<rfcavity>(*it)->getHarmonicNumber();
        } else if ( std::strcmp((*it)->Type(), "thinrfcavity") == 0) {
            sstream << ", Freq="
                    << boost::static_pointer_cast<thinrfcavity>(*it)->getRadialFrequency()/(2.0*mconstants::pi) <<
                    ", Harmon="
                    << boost::static_pointer_cast<thinrfcavity>(*it)->getHarmonicNumber();
        }
        sstream << std::endl;
    }
    return sstream.str();
}

void
print_full_chef_beamline(BmlPtr beamline_sptr)
{
    std::cout << full_chef_beamline_as_string(beamline_sptr) << std::endl;
    std::cout.flush();
}

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle)
{
    double mass = reference_particle.get_mass();
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
                    if (floating_point_equal(mass, PH_NORM_mmu,
                            mass_tolerance)) {
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
    Reference_particle reference_particle(
            static_cast<int >(chef_particle.Charge()), four_momentum);
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

void
ensure_jet_environment(int map_order)
{
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    } else if (Jet__environment::getLastEnv()->maxWeight() != map_order) {
        JetParticle::createStandardEnvironments(map_order);
    }
}

JetParticle
reference_particle_to_chef_jet_particle(
        Reference_particle const& reference_particle, int map_order)
{
    ensure_jet_environment(map_order);
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

Particle
get_closed_orbit_particle(Particle util_part, BmlPtr beamline_sptr, double dpop)
{
    BmlPtr throwaway_sptr(beamline_sptr->Clone());
    throwaway_sptr->setLineMode(beamline::ring);

    Jet__environment_ptr storedEnv = Jet__environment::getLastEnv();
    JetC__environment_ptr storedEnvC = JetC__environment::getLastEnv();
    // just need a basic first order environment
    JetParticle::createStandardEnvironments(1);

    ClosedOrbitSage closed_orbit_sage(throwaway_sptr);
    Particle probe(util_part);
    probe.set_ndp(dpop);
    JetParticle jetprobe(probe);
    closed_orbit_sage.findClosedOrbit(jetprobe);
    Particle closed_orbit_particle(jetprobe);
    // restore environment
    Jet__environment::setLastEnv(storedEnv);
    JetC__environment::setLastEnv(storedEnvC);
    return closed_orbit_particle;
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

