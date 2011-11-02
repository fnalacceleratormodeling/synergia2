#include "lattice_simulator.h"
#include <physics_toolkit/BeamlineContext.h>
#include <beamline/beamline_elements.h>
#include <basic_toolkit/PhysicsConstants.h>
#include <stdexcept>
#include <cstring>

Lattice_functions::Lattice_functions() :
    alpha_x(0.0), alpha_y(0.0), beta_x(0.0), beta_y(0.0), psi_x(0.0),
            psi_y(0.0), D_x(0.0), D_y(0.0), Dprime_x(0.0), Dprime_y(0.0)
{
}

Lattice_functions::Lattice_functions(LattFuncSage::lattFunc const& latt_func) :
    alpha_x(latt_func.alpha.hor), alpha_y(latt_func.alpha.ver),
            beta_x(latt_func.beta.hor), beta_y(latt_func.beta.ver),
            psi_x(latt_func.psi.hor), psi_y(latt_func.psi.ver),
            D_x(latt_func.dispersion.hor), D_y(latt_func.dispersion.ver),
            Dprime_x(latt_func.dPrime.hor), Dprime_y(latt_func.dPrime.ver)
{
}

void
Lattice_simulator::construct_extractor_map()
{
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    extractor_map_sptr->set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(
            chef_propagate_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_propagate_operation_extractor(chef_lattice_sptr,
                            map_order)));
    extractor_map_sptr->set_extractor(
            chef_map_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_map_operation_extractor(chef_lattice_sptr,
                            map_order)));
}

void
Lattice_simulator::get_tunes()
{
    if (!have_tunes) {
        BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
        BeamlineContext beamline_context(
                reference_particle_to_chef_particle(
                        lattice_sptr->get_reference_particle()), beamline_sptr);
        horizontal_tune = beamline_context.getHorizontalFracTune();
        vertical_tune = beamline_context.getVerticalFracTune();
        have_tunes = true;
    }
}

Lattice_simulator::Lattice_simulator(Lattice_sptr lattice_sptr, int map_order) :
    lattice_sptr(lattice_sptr),
            chef_lattice_sptr(new Chef_lattice(lattice_sptr)),
            extractor_map_sptr(new Operation_extractor_map),
            map_order(map_order), have_element_lattice_functions(false),
            have_slice_lattice_functions(false), have_tunes(false)
{
    construct_extractor_map();
}

void
Lattice_simulator::construct_sliced_chef_beamline(
        Lattice_element_slices const& slices)
{
    chef_lattice_sptr->construct_sliced_beamline(slices);
}

int
Lattice_simulator::get_map_order() const
{
    return map_order;
}

Operation_extractor_map_sptr
Lattice_simulator::get_operation_extractor_map_sptr()
{
    return extractor_map_sptr;
}

Lattice_sptr
Lattice_simulator::get_lattice_sptr()
{
    return lattice_sptr;
}

Chef_lattice_sptr
Lattice_simulator::get_chef_lattice_sptr()
{
    return chef_lattice_sptr;
}

void
Lattice_simulator::calculate_element_lattice_functions()
{
    if (!have_element_lattice_functions) {
        BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
        BeamlineContext beamline_context(
                reference_particle_to_chef_particle(
                        lattice_sptr->get_reference_particle()), beamline_sptr);
        std::vector<LattFuncSage::lattFunc > latt_func(
                beamline_context.getTwissArray());
        beamline::const_iterator it = beamline_sptr->begin();
        for (int i = 0; i < latt_func.size(); ++i) {
            ElmPtr chef_element(*it);
            if (chef_element->Name() != Chef_lattice::internal_marker_name) {
                Lattice_element lattice_element(
                        chef_lattice_sptr->get_lattice_element(chef_element));
                lattice_functions_element_map[&lattice_element] = latt_func.at(
                        i);
                ++it;
            }
        }
        have_element_lattice_functions = true;
    }
}

void
Lattice_simulator::calculate_slice_lattice_functions()
{
    if (!have_slice_lattice_functions) {
        BmlPtr beamline_sptr(chef_lattice_sptr->get_sliced_beamline_sptr());
        BeamlineContext beamline_context(
                reference_particle_to_chef_particle(
                        lattice_sptr->get_reference_particle()), beamline_sptr);
        std::vector<LattFuncSage::lattFunc > latt_func(
                beamline_context.getTwissArray());
        beamline::const_iterator it = beamline_sptr->begin();
        for (int i = 0; i < latt_func.size(); ++i) {
            ElmPtr chef_element(*it);
            if (chef_element->Name() != Chef_lattice::internal_marker_name) {
                Lattice_element_slice lattice_element_slice(
                        chef_lattice_sptr->get_lattice_element_slice(
                                chef_element));
                lattice_functions_slice_map[&lattice_element_slice]
                        = latt_func.at(i);
                ++it;
            }
        }
        have_slice_lattice_functions = true;
    }
}

Lattice_functions const&
Lattice_simulator::get_lattice_functions(Lattice_element & lattice_element)
{
    calculate_element_lattice_functions();
    return lattice_functions_element_map[&lattice_element];
}

Lattice_functions const&
Lattice_simulator::get_lattice_functions(
        Lattice_element_slice & lattice_element_slice)
{
    calculate_slice_lattice_functions();
    return lattice_functions_slice_map[&lattice_element_slice];
}

double
Lattice_simulator::get_horizontal_tune()
{
    get_tunes();
    return horizontal_tune;
}

double
Lattice_simulator::get_vertical_tune()
{
    get_tunes();
    return vertical_tune;
}

// set_chef_correctors is a local function
void
set_chef_correctors(Lattice_elements & correctors, Chef_lattice & chef_lattice,
        BeamlineContext & beamline_context, bool horizontal)
{
    for (Lattice_elements::iterator le_it = correctors.begin(); le_it
            != correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin(); ce_it
                != chef_elements.end(); ++ce_it) {
            if (std::strcmp((*ce_it)->Type(), "quadrupole") == 0) {
                if (horizontal) {
                    beamline_context.addHTuneCorrector(
                            boost::dynamic_pointer_cast<quadrupole >(*ce_it));
                } else {
                    beamline_context.addVTuneCorrector(
                            boost::dynamic_pointer_cast<quadrupole >(*ce_it));
                }
            } else if (std::strcmp((*ce_it)->Type(), "thinQuad") == 0) {
                if (horizontal) {
                    beamline_context.addHTuneCorrector(
                            boost::dynamic_pointer_cast<thinQuad >(*ce_it));
                } else {
                    beamline_context .addVTuneCorrector(
                            boost::dynamic_pointer_cast<thinQuad >(*ce_it));
                }
            } else {
                std::string message(
                        "Lattice_simulator::adjust_tunes: Lattice_element ");
                message += (*le_it)->get_name();
                message += " of type ";
                message += (*le_it)->get_type();
                message += " cannot be used as a corrector because it has a";
                message += " chef element of type ";
                message += (*ce_it)->Type();
                throw std::runtime_error(message.c_str());
            }
        }
    }
}

// extract_quad_strengths is a local function
void
extract_quad_strengths(Lattice_elements & correctors,
        Chef_lattice & chef_lattice)
{
    for (Lattice_elements::iterator le_it = correctors.begin(); le_it
            != correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin(); ce_it
                != chef_elements.end(); ++ce_it) {
            double k1 = (*ce_it)->Strength() / chef_lattice.get_brho();
            (*le_it)->set_double_attribute("k1", k1);
        }
    }
}

void
Lattice_simulator::adjust_tunes(double horizontal_tune, double vertical_tune,
        Lattice_elements & horizontal_correctors,
        Lattice_elements & vertical_correctors)
{
    BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
    BeamlineContext beamline_context(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()), beamline_sptr);
    set_chef_correctors(horizontal_correctors, *chef_lattice_sptr,
            beamline_context, true);
    set_chef_correctors(vertical_correctors, *chef_lattice_sptr,
            beamline_context, false);
    int status = beamline_context.changeTunesTo(horizontal_tune, vertical_tune);
    if (status == BeamlineContext::NO_TUNE_ADJUSTER) {
        throw std::runtime_error(
                "Lattice_simulator::adjust_tunes: no corrector elements found");
    } else if (status != BeamlineContext::OKAY) {
        throw std::runtime_error(
                "Lattice_simulator::adjust_tunes: failed with unknown status");
    }
    extract_quad_strengths(horizontal_correctors, *chef_lattice_sptr);
    extract_quad_strengths(vertical_correctors, *chef_lattice_sptr);
}

Lattice_simulator::~Lattice_simulator()
{
}

