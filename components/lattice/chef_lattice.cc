#include "chef_lattice.h"
#include "utils/floating_point.h"
#include <beamline/beamline_elements.h>
#include <basic_toolkit/PhysicsConstants.h>
#include <physics_toolkit/DriftConverter.h>
#include <beamline/RefRegVisitor.h>

#include <stdexcept>

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle)
{
    // n.b. We don't deal with negatively charged particles!!!!
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

beamline
Chef_lattice::construct_raw_lattice(Lattice_element_to_chef_fn_map const& map)
{
    beamline raw_beamline;
    for (Lattice_elements::const_iterator latt_it =
            lattice_ptr->get_elements().begin(); latt_it
            != lattice_ptr->get_elements().end(); ++latt_it) {
        Lattice_element_to_chef_fn_map::const_iterator map_it = map.find(
                latt_it->get_type());
        if (map_it == map.end()) {
            throw(runtime_error("Chef_lattice: " + latt_it->get_type()
                    + " not handled"));
        } else {
            Chef_element_list celm_list = map_it->second(*latt_it, brho);
            for (Chef_element_list::const_iterator cel_it = celm_list.begin(); cel_it
                    != celm_list.end(); ++cel_it) {
                raw_beamline.append(*cel_it);
            }
        }
    }
    return raw_beamline;
}

void
Chef_lattice::polish_lattice(beamline const& raw_beamline)
{
    DriftConverter drift_converter;
    beamline_ptr = drift_converter.convert(raw_beamline);
    Particle testpart(reference_particle_to_chef_particle(
            lattice_ptr->get_reference_particle()));
    RefRegVisitor registrar(testpart);
    beamline_ptr->accept(registrar);
}

void
Chef_lattice::construct(Lattice_element_to_chef_fn_map const& map)
{
    if (!lattice_ptr->has_reference_particle()) {
        throw(std::runtime_error(
                "Chef_lattice: requires a reference particle in Lattice"));
    }
    brho = lattice_ptr->get_reference_particle().get_momentum()
            / PH_CNV_brho_to_p;

    polish_lattice(construct_raw_lattice(map));
}

Chef_lattice::Chef_lattice(Lattice & lattice) :
    lattice_ptr(&lattice), beamline_ptr()
{
    construct(get_standard_lattice_element_to_chef_fn_map());
}

Chef_lattice::Chef_lattice(Lattice & lattice,
        Lattice_element_to_chef_fn_map const& map) :
    lattice_ptr(&lattice), beamline_ptr()
{
    construct(map);
}

BmlPtr
Chef_lattice::get_beamline_ptr()
{
    return beamline_ptr;
}

Chef_lattice::~Chef_lattice()
{

}

Lattice_element_to_chef_fn_map
get_standard_lattice_element_to_chef_fn_map()
{
    Lattice_element_to_chef_fn_map map;
    map["marker"] = lattice_element_to_chef_marker;
    map["drift"] = lattice_element_to_chef_drift;
    //    map["sbend"] = lattice_element_to_chef_sbend;
    //    map["rbend"] = lattice_element_to_chef_rbend;
    map["quadrupole"] = lattice_element_to_chef_quadrupole;
    //    map["sextupole"] = lattice_element_to_chef_sextupole;
    //    map["octupole"] = lattice_element_to_chef_octupole;
    //    map["multipole"] = lattice_element_to_chef_multipole;
    //    map["solenoid"] = lattice_element_to_chef_solenoid;
    //    map["hkicker"] = lattice_element_to_chef_hkicker;
    //    map["vkicker"] = lattice_element_to_chef_vkicker;
    //    map["kicker"] = lattice_element_to_chef_kicker;
    //    map["rfcavity"] = lattice_element_to_chef_rfcavity;
    //    map["elseparator"] = lattice_element_to_chef_elseperator;
    //    map["hmonitor"] = lattice_element_to_chef_hmonitor;
    //    map["vmonitor"] = lattice_element_to_chef_vmonitor;
    //    map["monitor"] = lattice_element_to_chef_monitor;
    //    map["instrument"] = lattice_element_to_chef_instrument;
    //    map["ecollimator"] = lattice_element_to_chef_ecollimator;
    //    map["rcollimator"] = lattice_element_to_chef_rcollimator;

    return map;
}

Chef_element_list
lattice_element_to_chef_marker(Lattice_element const& lattice_element,
        double brho)
{
    Chef_element_list retval;
    ElmPtr elm(new marker(lattice_element.get_name().c_str()));
    retval.push_back(elm);
    return retval;
}

Chef_element_list
lattice_element_to_chef_drift(Lattice_element const& lattice_element,
        double brho)
{
    Chef_element_list retval;
    ElmPtr elm(new drift(lattice_element.get_name().c_str(),
            lattice_element.get_length()));
    retval.push_back(elm);
    return retval;
}

Chef_element_list
lattice_element_to_chef_quadrupole(Lattice_element const& lattice_element,
        double brho)
{
    // tilt can have a string value of ""
    if (lattice_element.has_string_attribute("tilt")) {
        throw(runtime_error(
                "lattice_element_to_chef_quadrupole: tilt element not handled"));
    }
    if (lattice_element.has_double_attribute("tilt")) {
        if (lattice_element.get_double_attribute("tilt") != 0.0) {
            throw(runtime_error(
                    "lattice_element_to_chef_quadrupole: non-zero tilt element not handled"));
        }
    }
    Chef_element_list retval;

    double length = lattice_element.get_length();
    bmlnElmnt* bmln_elmnt;
    if (length == 0.0) {
        bmln_elmnt = new quadrupole(lattice_element.get_name().c_str(), length,
                brho * lattice_element.get_double_attribute("k1"));
    } else {
        bmln_elmnt = new thinQuad(lattice_element.get_name().c_str(), brho
                * lattice_element.get_double_attribute("k1"));
    }
    ElmPtr elm(bmln_elmnt);
    retval.push_back(elm);
    return retval;
}
