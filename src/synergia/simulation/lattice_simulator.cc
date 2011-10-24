#include "lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include <physics_toolkit/BeamlineContext.h>
#include <beamline/beamline_elements.h>
#include <basic_toolkit/PhysicsConstants.h>
#include <stdexcept>
#include <cstring>

Lattice_functions::Lattice_functions() :
    alpha_x(0.0), alpha_y(0.0), beta_x(0.0), beta_y(0.0), psi_x(0.0),
            psi_y(0.0), D_x(0.0), D_y(0.0), Dprime_x(0.0), Dprime_y(0.0),
            arc_length(0.0)
{
}

Lattice_functions::Lattice_functions(LattFuncSage::lattFunc const& latt_func) :
    alpha_x(latt_func.alpha.hor), alpha_y(latt_func.alpha.ver),
            beta_x(latt_func.beta.hor), beta_y(latt_func.beta.ver),
            psi_x(latt_func.psi.hor), psi_y(latt_func.psi.ver),
            D_x(latt_func.dispersion.hor), D_y(latt_func.dispersion.ver),
            Dprime_x(latt_func.dPrime.hor), Dprime_y(latt_func.dPrime.ver),
            arc_length(latt_func.arcLength)
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
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
    if (!have_tunes) {
        BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
        BeamlineContext beamline_context(
                reference_particle_to_chef_particle(
                        lattice_sptr->get_reference_particle()), beamline_sptr);
        horizontal_tune = beamline_context.getHorizontalFracTune();
        vertical_tune = beamline_context.getVerticalFracTune();
        have_tunes = true;
    }

void
Lattice_simulator::construct_aperture_extractor_map()
{
    typedef Generic_aperture_extractor<Circular_aperture_operation >
            Circular_extractor;
    aperture_extractor_map_sptr->set_extractor(
            Circular_aperture_operation::attribute_name,
            boost::shared_ptr<Circular_extractor >(new Circular_extractor()));

    typedef Generic_aperture_extractor<Elliptical_aperture_operation >
            Elliptical_extractor;
    aperture_extractor_map_sptr->set_extractor(
            Elliptical_aperture_operation::attribute_name,
            boost::shared_ptr<Elliptical_extractor >(new Elliptical_extractor()));

    typedef Generic_aperture_extractor<Rectangular_aperture_operation >
            Rectangular_extractor;
    aperture_extractor_map_sptr->set_extractor(
            Rectangular_aperture_operation::attribute_name,
            boost::shared_ptr<Rectangular_extractor >(
                    new Rectangular_extractor()));
}

Lattice_simulator::Lattice_simulator(Lattice_sptr lattice_sptr, int map_order) :
    lattice_sptr(lattice_sptr),
            chef_lattice_sptr(new Chef_lattice(lattice_sptr)),
            extractor_map_sptr(new Operation_extractor_map),
            have_element_lattice_functions(false),
            have_slice_lattice_functions(false), have_tunes(false)
            aperture_extractor_map_sptr(new Aperture_operation_extractor_map),
            map_order(map_order), have_slices(false),
{
    construct_extractor_map();
    construct_aperture_extractor_map();
    set_bucket_length();
}

void
Lattice_simulator::set_slices(Lattice_element_slices const& slices)
{
    this->slices = slices;
    have_slices = true;
    construct_sliced_chef_beamline();
}

void
Lattice_simulator::construct_sliced_chef_beamline()
{
    if (!have_slices) {
        throw std::runtime_error(
                "Lattice_simulator::construct_sliced_chef_beamline called before set_slices");
    }
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

Aperture_operation_extractor_map_sptr
Lattice_simulator::get_aperture_operation_extractor_map_sptr()
{
    return aperture_extractor_map_sptr;
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
Lattice_simulator::set_bucket_length()
{
    double freq(0.), freq2(0.);
    int isw = 0;
    double eps = 1e-6;
    for (Lattice_elements::const_iterator it =
            this->lattice_sptr->get_elements().begin(); it
            != this->lattice_sptr->get_elements().end(); ++it) {

        if ((*it)->has_double_attribute("freq")) {
            freq = (*it)->get_double_attribute("freq");
            if ((isw == 1) && (fabs(freq - freq2) > eps)) {
                throw std::runtime_error(
                        "set_bucket_length: rf elements with different frequencies found!!");
            }
            freq2 = freq;
            isw = 1;
        }
        if (isw == 1) {
            double
                    beta =
                            this->get_lattice_sptr()->get_reference_particle().get_beta();
            this->bucket_length = pconstants::c * beta / freq;
        } else {
            this->bucket_length = 0.0;
        }
    }
}

double
Lattice_simulator::get_bucket_length()
{
    return this->bucket_length;
}

int
Lattice_simulator::get_number_buckets()
{
    double eps = 1e-5;
    int number_buckets;
    double bl = get_bucket_length();
    double ol = this->get_lattice_sptr()->get_length();
    bl > eps ? number_buckets = int(ol / bl) : number_buckets = 1;
    return number_buckets;
}

void
Lattice_simulator::update()
{
    chef_lattice_sptr = Chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    construct_extractor_map();
    construct_aperture_extractor_map();
    if (have_slices) {
        construct_sliced_chef_beamline();
    }
    have_element_lattice_functions = false;
    lattice_functions_element_map.clear();
    have_slice_lattice_functions = false;
    lattice_functions_slice_map.clear();
    have_tunes = false;
}

void
Lattice_simulator::calculate_element_lattice_functions()
{
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
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
            if (std::strcmp(chef_element->Name().c_str(),
                    Chef_lattice::internal_marker_name)) {
                Lattice_element const
                        * lattice_element_ptr(
                                &(chef_lattice_sptr->get_lattice_element(
                                        chef_element)));
                lattice_functions_element_map[lattice_element_ptr]
                        = Lattice_functions(latt_func.at(i));
            }
            ++it;
        }
        have_element_lattice_functions = true;
    }
}

void
Lattice_simulator::calculate_slice_lattice_functions()
{
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
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
            if (std::strcmp(chef_element->Name().c_str(),
                    Chef_lattice::internal_marker_name)) {
                Lattice_element_slice const * lattice_element_slice_ptr(
                        &(chef_lattice_sptr->get_lattice_element_slice(
                                chef_element)));
                lattice_functions_slice_map[lattice_element_slice_ptr]
                        = latt_func.at(i);
            }
            ++it;
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
set_chef_correctors(Lattice_elements const& correctors,
        Chef_lattice & chef_lattice, BeamlineContext & beamline_context,
        bool horizontal)
{
    for (Lattice_elements::const_iterator le_it = correctors.begin(); le_it
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
extract_quad_strengths(Lattice_elements const& correctors,
        Chef_lattice & chef_lattice)
{
    for (Lattice_elements::const_iterator le_it = correctors.begin(); le_it
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
        Lattice_elements const& horizontal_correctors,
        Lattice_elements const& vertical_correctors, double tolerance)
{
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
    BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
    BeamlineContext beamline_context(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()), beamline_sptr);
    set_chef_correctors(horizontal_correctors, *chef_lattice_sptr,
            beamline_context, true);
    set_chef_correctors(vertical_correctors, *chef_lattice_sptr,
            beamline_context, false);
    double nu_h = beamline_context.getHorizontalFracTune();
    double nu_v = beamline_context.getVerticalFracTune();
    double horizontal_final_tune = horizontal_tune - double(int(horizontal_tune));
    double vertical_final_tune = vertical_tune - double(int(vertical_tune));
    int rank = Commxx().get_rank();
    if (rank == 0) {
        std::cout << "        Initial tunes: horizontal    : " << nu_h
                << std::endl;;
        std::cout << "                       vertical      : " << nu_v
                << std::endl;
        std::cout << "        Final tunes: horizontal      : "
                << horizontal_final_tune << std::endl;
        std::cout << "                     vertical        : "
                << vertical_final_tune << std::endl;
    }
#if 1
    double adjuster_tune_step = 0.105;
    double const dtune_h =  (horizontal_final_tune - nu_h) / double(1 + int((
            std::abs(horizontal_final_tune - nu_h)) / adjuster_tune_step));

    double horizontal_target_tune = nu_h + dtune_h;
    double vertical_target_tune = vertical_final_tune;
    int step = 1;
    while ((10.0 * tolerance) < std::abs(horizontal_final_tune - nu_h)) {
        if (rank == 0) std::cout << "\n        Step " << step << std::endl;
        while (tolerance < std::abs( horizontal_target_tune - nu_h)) {
            int status = beamline_context.changeTunesTo(horizontal_target_tune,
                     vertical_target_tune);
            if (status == BeamlineContext::NO_TUNE_ADJUSTER) {
                throw std::runtime_error(
                        "Lattice_simulator::adjust_tunes: no corrector elements found");
            } else if (status != BeamlineContext::OKAY) {
                throw std::runtime_error(
                        "Lattice_simulator::adjust_tunes: failed with unknown status");
        }
            nu_h = beamline_context.getHorizontalFracTune();
            nu_v = beamline_context.getVerticalFracTune();
            if (rank == 0) {
                std::cout << "        Tunes: horizontal: " << nu_h << " = "
                        << horizontal_target_tune << " + "
                        << nu_h - horizontal_target_tune << std::endl;
                std::cout << "               vertical  : " << nu_v << std::endl;
            }
        }
        have_tunes = false;
        horizontal_target_tune += dtune_h;
        ++step;
    }
#endif
#if 0
    const int max_iterations = 20;
    int iteration = 0;
    while ((std::abs(get_horizontal_tune() - horizontal_tune) > tolerance)
            && (std::abs(get_vertical_tune() - vertical_tune) > tolerance)) {
        ++iteration;
        if (iteration > max_iterations) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_tunes: maximum number of iterations exceeded");
        }
        int status = beamline_context.changeTunesTo(horizontal_tune,
                vertical_tune);
        if (status == BeamlineContext::NO_TUNE_ADJUSTER) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_tunes: no corrector elements found");
        } else if (status != BeamlineContext::OKAY) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_tunes: failed with unknown status");
        }
        have_tunes = false;
    }
#endif
    extract_quad_strengths(horizontal_correctors, *chef_lattice_sptr);
    extract_quad_strengths(vertical_correctors, *chef_lattice_sptr);
    update();
}

Lattice_simulator::~Lattice_simulator()
{
}

