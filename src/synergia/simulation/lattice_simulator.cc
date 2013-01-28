#include "lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"

#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wsequence-point"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#include <beamline/beamline_elements.h>
#include <physics_toolkit/Sage.h>
#include <basic_toolkit/PhysicsConstants.h>
#include <beamline/RefRegVisitor.h>
#if __GNUC__ > 4 && __GNUC_MINOR__ > 5
#pragma GCC diagnostic pop
#endif

#include <stdexcept>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>

Lattice_functions::Lattice_functions() :
        alpha_x(0.0), alpha_y(0.0), beta_x(0.0), beta_y(0.0), psi_x(0.0), psi_y(
                0.0), D_x(0.0), D_y(0.0), Dprime_x(0.0), Dprime_y(0.0), arc_length(
                0.0)
{
}

Lattice_functions::Lattice_functions(LattFuncSage::lattFunc const& latt_func) :
        alpha_x(latt_func.alpha.hor), alpha_y(latt_func.alpha.ver), beta_x(
                latt_func.beta.hor), beta_y(latt_func.beta.ver), psi_x(
                latt_func.psi.hor), psi_y(latt_func.psi.ver), D_x(
                latt_func.dispersion.hor), D_y(latt_func.dispersion.ver), Dprime_x(
                latt_func.dPrime.hor), Dprime_y(latt_func.dPrime.ver), arc_length(
                latt_func.arcLength)
{
}

template<class Archive>
    void
    Lattice_functions::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(alpha_x);
        ar & BOOST_SERIALIZATION_NVP(alpha_y);
        ar & BOOST_SERIALIZATION_NVP(beta_x);
        ar & BOOST_SERIALIZATION_NVP(beta_y);
        ar & BOOST_SERIALIZATION_NVP(psi_x);
        ar & BOOST_SERIALIZATION_NVP(psi_y);
        ar & BOOST_SERIALIZATION_NVP(D_x);
        ar & BOOST_SERIALIZATION_NVP(D_y);
        ar & BOOST_SERIALIZATION_NVP(Dprime_x);
        ar & BOOST_SERIALIZATION_NVP(Dprime_y);
        ar & BOOST_SERIALIZATION_NVP(arc_length);
    }

template
void
Lattice_functions::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lattice_functions::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lattice_functions::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lattice_functions::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

void
Lattice_simulator::construct_extractor_map()
{
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    extractor_map_sptr->set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_propagate_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_propagate_operation_extractor(chef_lattice_sptr,
                            map_order)));
    extractor_map_sptr->set_extractor(chef_map_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_map_operation_extractor(chef_lattice_sptr,
                            map_order)));
}

void
Lattice_simulator::calculate_beamline_context()
{
    ensure_jet_environment(map_order);
    BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
    beamline_context_sptr = BmlContextPtr(
            new BeamlineContext(
                    reference_particle_to_chef_particle(
                            lattice_sptr->get_reference_particle()),
                    beamline_sptr));
    if (!Sage::isRing(beamline_sptr)) {
        beamline_context_sptr->handleAsRing();
    }
}

void
Lattice_simulator::calculate_sliced_beamline_context()
{
    if (!have_slices) {
        throw std::runtime_error(
                "Lattice_simulator::calculate_sliced_beamline_context called before set_slices");
    }
    ensure_jet_environment(map_order);
    BmlPtr beamline_sptr(chef_lattice_sptr->get_sliced_beamline_sptr());
    sliced_beamline_context_sptr = BmlContextPtr(
            new BeamlineContext(
                    reference_particle_to_chef_particle(
                            lattice_sptr->get_reference_particle()),
                    beamline_sptr));
    if (!Sage::isRing(beamline_sptr)) {
        sliced_beamline_context_sptr->handleAsRing();
    }
}

BmlContextPtr
Lattice_simulator::get_beamline_context()
{
    if (!have_beamline_context) {
        calculate_beamline_context();
        have_beamline_context = true;
    }
    return beamline_context_sptr;
}

BmlContextPtr
Lattice_simulator::get_sliced_beamline_context()
{
    if (!have_sliced_beamline_context) {
        calculate_sliced_beamline_context();
        have_sliced_beamline_context = true;
    }
    return sliced_beamline_context_sptr;
}

bool
Lattice_simulator::is_ring()
{
    get_beamline_context();
    return (beamline_context_sptr->isRing());
}

void
Lattice_simulator::get_tunes()
{
    if (!have_tunes) {
        get_beamline_context();
        horizontal_tune = beamline_context_sptr->getHorizontalFracTune();
        vertical_tune = beamline_context_sptr->getVerticalFracTune();
        have_tunes = true;
    }
}

void
Lattice_simulator::construct_aperture_extractor_map()
{
    aperture_extractor_map_sptr->set_extractor(
            Circular_aperture_operation::attribute_name,
            boost::shared_ptr<Circular_extractor >(new Circular_extractor()));

    aperture_extractor_map_sptr->set_extractor("default",
            boost::shared_ptr<Circular_extractor >(new Circular_extractor()));

    aperture_extractor_map_sptr->set_extractor(
            Elliptical_aperture_operation::attribute_name,
            boost::shared_ptr<Elliptical_extractor >(
                    new Elliptical_extractor()));

    aperture_extractor_map_sptr->set_extractor(
            Rectangular_aperture_operation::attribute_name,
            boost::shared_ptr<Rectangular_extractor >(
                    new Rectangular_extractor()));

    aperture_extractor_map_sptr->set_extractor(
            Polygon_aperture_operation::attribute_name,
            boost::shared_ptr<Polygon_extractor >(new Polygon_extractor()));

    aperture_extractor_map_sptr->set_extractor(
            Wire_elliptical_aperture_operation::attribute_name,
            boost::shared_ptr<Wire_elliptical_extractor >(
                    new Wire_elliptical_extractor()));

    aperture_extractor_map_sptr->set_extractor(
            Lambertson_aperture_operation::attribute_name,
            boost::shared_ptr<Lambertson_extractor >(
                    new Lambertson_extractor()));
}

Lattice_simulator::Lattice_simulator(Lattice_sptr lattice_sptr, int map_order) :
        lattice_sptr(lattice_sptr), have_slices(false), chef_lattice_sptr(
                new Chef_lattice(lattice_sptr)), extractor_map_sptr(
                new Operation_extractor_map), aperture_extractor_map_sptr(
                new Aperture_operation_extractor_map), have_beamline_context(
                false), have_sliced_beamline_context(false), map_order(
                map_order), have_element_lattice_functions(false), have_slice_lattice_functions(
                false), horizontal_tune(0.0), vertical_tune(0.0), have_tunes(
                false), horizontal_chromaticity(0.0), vertical_chromaticity(
                0.0), have_chromaticities(false), linear_one_turn_map(
                boost::extents[6][6])
{
    construct_extractor_map();
    construct_aperture_extractor_map();
    set_bucket_length();
}

Lattice_simulator::Lattice_simulator()
{
}

void
Lattice_simulator::set_slices(Lattice_element_slices const& slices)
{
    this->slices = slices;
    have_slices = true;
    construct_sliced_chef_beamline();
}

Lattice_element_slices const&
Lattice_simulator::get_slices() const
{
    if (!have_slices) {
        throw std::runtime_error(
                "Lattice_simulator::get_slices called before set_slices");
    }
    return slices;
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
            this->lattice_sptr->get_elements().begin();
            it != this->lattice_sptr->get_elements().end(); ++it) {

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
            double beta =
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
    have_beamline_context = false;
    have_chromaticities = false;
    normal_form_sage_sptr.reset();
}

void
Lattice_simulator::calculate_element_lattice_functions()
{
    if (!have_element_lattice_functions) {
        BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
        std::vector<LattFuncSage::lattFunc > latt_func(
                get_beamline_context()->getTwissArray());
        beamline::const_iterator it = beamline_sptr->begin();
        for (unsigned int i = 0; i < latt_func.size(); ++i) {
            ElmPtr chef_element(*it);
            if (std::strcmp(chef_element->Name().c_str(),
                    Chef_lattice::internal_marker_name)) {
                Lattice_element * lattice_element_ptr(
                        &(chef_lattice_sptr->get_lattice_element(chef_element)));
                lattice_functions_element_map[lattice_element_ptr] =
                        Lattice_functions(latt_func.at(i));
            }
            ++it;
        }
        have_element_lattice_functions = true;
    }
}

void
Lattice_simulator::calculate_slice_lattice_functions()
{
    if (!have_slice_lattice_functions) {
        BmlPtr beamline_sptr(chef_lattice_sptr->get_sliced_beamline_sptr());
        std::vector<LattFuncSage::lattFunc > latt_func(
                get_sliced_beamline_context()->getTwissArray());
        beamline::const_iterator it = beamline_sptr->begin();
        for (unsigned int i = 0; i < latt_func.size(); ++i) {
            ElmPtr chef_element(*it);
            if (std::strcmp(chef_element->Name().c_str(),
                    Chef_lattice::internal_marker_name)) {
                Lattice_element_slice * lattice_element_slice_ptr(
                        &(chef_lattice_sptr->get_lattice_element_slice(
                                chef_element)));
                lattice_functions_slice_map[lattice_element_slice_ptr] =
                        latt_func.at(i);
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

// calculate the normal form sage object, leave a shared pointer to this object
// in the lattice_simulator class member normal_form_sage_sptr
void
Lattice_simulator::calculate_normal_form()
{
    get_beamline_context();
    Mapping one_turn_map = beamline_context_sptr->getOneTurnMap();
    normal_form_sage_sptr = Normal_form_sage_sptr(
            new normalFormSage(one_turn_map,
                    reference_particle_to_chef_jet_particle(
                            lattice_sptr->get_reference_particle(), map_order),
                    map_order));
}

// return the normal_form_sage_sptr if it exists, otherwise calculate
// it first and then return it.
Normal_form_sage_sptr
Lattice_simulator::get_normal_form_sptr()
{
    if (normal_form_sage_sptr) {
        return normal_form_sage_sptr;
    } else {
        calculate_normal_form();
        return normal_form_sage_sptr;
    }
}

// checkLinearNormalForm
// check the linear part of the normal form calculation
bool
Lattice_simulator::check_linear_normal_form()
{
    Normal_form_sage_sptr nf_sptr(get_normal_form_sptr());
    return nf_sptr->checkLinearNormalForm();
}

// converts a MArray2d of particle coordinates in synergia ordering into
// an MArray2d of complex normal form coordinates stored as a0.real, a0.imag,
//  a1.real,a1.imag, a2.real, a2.imag.

void
Lattice_simulator::convert_human_to_normal(MArray2d_ref coords)
{
    Normal_form_sage_sptr nf_sptr(get_normal_form_sptr());

#if 0
    std::cout << "chef->synergia indices:" << std::endl;
    for (int j=0; j<6; ++j) {
        std::cout << j << " ->  " << get_synergia_index(j) << std::endl;
    }
    std::cout << std::endl;
    std::cout << "synergia->chef indices:" << std::endl;
    for (int j=0; j<6; ++j) {
        std::cout << j << " ->  " << get_chef_index(j) << std::endl;
    }
#endif
    const MArray2d::size_type *coords_shape = coords.shape();
    const MArray2d::index *coords_bases = coords.index_bases();

#if 0
    std::cout << "coords_shape: " << coords_shape[0] << ", " << coords_shape[1] << std::endl;
    std::cout << "coords_bases: " << coords_bases[0] << ", " << coords_bases[1] << std::endl;
#endif

    if ((coords_shape[1] != 7) || (coords_bases[1] != 0)) {
        throw std::runtime_error(
                "Lattice_simulator::convert_human_to_normal expected nx[0:7] array");
    }

    for (unsigned int i = coords_bases[0];
            i != coords_bases[0] + coords_shape[0]; ++i) {
        Vector w(6);
        VectorC a(6);

        for (int j = 0; j < 6; ++j) {
            w(get_chef_index(j)) = coords[i][j];
        }
#if 0
        std::cout << "human->normal human(chef): " << w << std::endl;
#endif
        nf_sptr->cnvDataToNormalForm(w, a);
#if 0
        std::cout << "human->normal normal(chef): " << a << std::endl;
#endif

        for (int j = 0; j < 3; ++j) {
            coords[i][2 * j] = a(j).real();
            coords[i][2 * j + 1] = a(j).imag();
        }
    }
}

// converts a MArray2d of complex normal form particle coordinates into
// human space coordinates in synergia order.
void
Lattice_simulator::convert_normal_to_human(MArray2d_ref coords)
{
    Normal_form_sage_sptr nf_sptr(get_normal_form_sptr());

    const MArray2d::size_type *coords_shape = coords.shape();
    const MArray2d::index *coords_bases = coords.index_bases();

    if ((coords_shape[1] != 7) || (coords_bases[1] != 0)) {
        throw std::runtime_error(
                "Lattice_simulator::convert_normal_to_human expected nx[0:7] array");
    }
    for (unsigned int i = coords_bases[0];
            i != coords_bases[0] + coords_shape[0]; ++i) {
        Vector w(6);
        VectorC a(6);

        for (int j = 0; j < 3; ++j) {
            a(j) = std::complex<double >(coords[i][2 * j],
                    coords[i][2 * j + 1]);
            a(j + 3) = std::conj(a(j));
        }

        // convert to human form in CHEF order
        nf_sptr->cnvDataFromNormalForm(a, w);

        // write back into synergia ordering
        for (int j = 0; j < 6; ++j) {
            coords[i][get_synergia_index(j)] = w(j);
        }
    }
}

// return a vector of the mean actions that will generate a stationary
// beam distribution having the specified standard deviations in each
// of the three planes.
std::vector<double >
Lattice_simulator::get_stationary_actions(const double stdx, const double stdy,
        const double std_cdt)
{
    Normal_form_sage_sptr nf_sptr(get_normal_form_sptr());
    // stationaryActions wants the second moments of the canonical variables
    // which are x,y,t.  Note that the input longitudinal variaable is std_cdt.
    // Convert that into stdt.
    double stdt = std_cdt / pconstants::c;

    std::vector<double > v(nf_sptr->stationaryActions(stdx, stdy, stdt));
    return v;
}

// returns the linear one turn map for the lattice and beam parameters
// for this lattice_simulator
Const_MArray2d_ref
Lattice_simulator::get_linear_one_turn_map()
{

    get_beamline_context();
    MatrixD lin_one_turn_map =
            beamline_context_sptr->getOneTurnMap().Jacobian();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            linear_one_turn_map[i][j] = lin_one_turn_map(get_chef_index(i),
                    get_chef_index(j));
        }
    }
    return linear_one_turn_map;
}

// this method should be used instead of the get_horizontal_tune() and get_vertical_tune()
std::pair<double, double >
Lattice_simulator::get_both_tunes()
{
    get_tunes();
    update(); // remake CHEF beamline to restore RF turned off by getHorizontalFracTune()
    return std::pair<double, double >(horizontal_tune, vertical_tune);
}

double
Lattice_simulator::get_horizontal_tune()
{
    get_tunes();
    update(); // remake CHEF beamline to restore RF turned off by getHorizontalFracTune()
    return horizontal_tune;
}

double
Lattice_simulator::get_vertical_tune()
{
    get_tunes();
    update(); // remake CHEF beamline to restore RF turned off by getHorizontalFracTune()
    return vertical_tune;
}

// set_chef_correctors is a local function
void
set_chef_correctors(Lattice_elements const& correctors,
        Chef_lattice & chef_lattice, BmlContextPtr & beamline_context_sptr,
        bool horizontal)
{
    for (Lattice_elements::const_iterator le_it = correctors.begin();
            le_it != correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin();
                ce_it != chef_elements.end(); ++ce_it) {
            if (std::strcmp((*ce_it)->Type(), "quadrupole") == 0) {
                if (horizontal) {
                    beamline_context_sptr->addHTuneCorrector(
                            boost::dynamic_pointer_cast<quadrupole >(*ce_it));
                } else {
                    beamline_context_sptr->addVTuneCorrector(
                            boost::dynamic_pointer_cast<quadrupole >(*ce_it));
                }
            } else if (std::strcmp((*ce_it)->Type(), "thinQuad") == 0) {
                if (horizontal) {
                    beamline_context_sptr->addHTuneCorrector(
                            boost::dynamic_pointer_cast<thinQuad >(*ce_it));
                } else {
                    beamline_context_sptr->addVTuneCorrector(
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
    for (Lattice_elements::const_iterator le_it = correctors.begin();
            le_it != correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin();
                ce_it != chef_elements.end(); ++ce_it) {
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
    BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
    get_beamline_context();
    set_chef_correctors(horizontal_correctors, *chef_lattice_sptr,
            beamline_context_sptr, true);
    set_chef_correctors(vertical_correctors, *chef_lattice_sptr,
            beamline_context_sptr, false);
    double nu_h = beamline_context_sptr->getHorizontalFracTune();
    double nu_v = beamline_context_sptr->getVerticalFracTune();
    double horizontal_final_tune = horizontal_tune
            - double(int(horizontal_tune));
    double vertical_final_tune = vertical_tune - double(int(vertical_tune));
    int rank = Commxx().get_rank();
    if (rank == 0) {
        std::cout << "        Initial tunes: horizontal    : " << nu_h
                << std::endl;
        ;
        std::cout << "                       vertical      : " << nu_v
                << std::endl;
        std::cout << "        Final tunes: horizontal      : "
                << horizontal_final_tune << std::endl;
        std::cout << "                     vertical        : "
                << vertical_final_tune << std::endl;
    }

    int step = 0;
    bool in_tolerance = sqrt(
            pow((nu_h - horizontal_final_tune), 2)
                    + pow((nu_v - vertical_final_tune), 2)) < tolerance;

    while (!in_tolerance && (step < 20)) {
        if (rank == 0) std::cout << "\n        Step " << step << std::endl;
        int status = beamline_context_sptr->changeTunesTo(horizontal_final_tune,
                vertical_final_tune);
        if (status == BeamlineContext::NO_TUNE_ADJUSTER) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_tunes: no corrector elements found");
        } else if (status != BeamlineContext::OKAY) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_tunes: failed with unknown status");
        }
        nu_h = beamline_context_sptr->getHorizontalFracTune();
        nu_v = beamline_context_sptr->getVerticalFracTune();
        if (rank == 0) {
            std::cout << "        Step tunes: horizontal: " << nu_h
                    << ", error " << horizontal_final_tune - nu_h << std::endl;
            std::cout << "                    vertical  : " << nu_v
                    << ", error " << vertical_final_tune - nu_v << std::endl;
        }
        in_tolerance = sqrt(
                pow((nu_h - horizontal_final_tune), 2)
                        + pow((nu_v - vertical_final_tune), 2)) < tolerance;
        ++step;
    }
    if (!in_tolerance) {
        std::cout << "Error, did not meet final tolerance within 20 steps"
                << std::endl;
    }
    extract_quad_strengths(horizontal_correctors, *chef_lattice_sptr);
    extract_quad_strengths(vertical_correctors, *chef_lattice_sptr);
    update();
}

void
Lattice_simulator::get_chromaticities()
{

    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
    if (!have_chromaticities) {
        double momentum(lattice_sptr->get_reference_particle().get_momentum());
        Proton probe;
        probe.SetReferenceMomentum(momentum);
        probe.setStateToZero();

        BmlPtr copy_beamline_sptr(
                chef_lattice_sptr->get_beamline_sptr()->Clone());
        copy_beamline_sptr->setEnergy(probe.ReferenceEnergy());

        BeamlineContext probecontext(probe, copy_beamline_sptr);
        probecontext.handleAsRing();

        double hcentral_tune, vcentral_tune;
        //  hcentral_tune = probecontext.getHorizontalFracTune();
        //  vcentral_tune = probecontext.getVerticalFracTune();

        hcentral_tune = probecontext.getHorizontalEigenTune();
        vcentral_tune = probecontext.getVerticalEigenTune();

        double chromat_H = 0.;
        double chromat_V = 0.;
        double dppcount = 0.;
        for (double dpp = -0.0005; dpp <= 0.00051; dpp += 0.0002) {
            Proton newprobe;
            newprobe.SetReferenceMomentum(momentum * (1.0 + dpp));
            newprobe.setStateToZero();
            copy_beamline_sptr->setEnergy(newprobe.ReferenceEnergy());
            BeamlineContext probecontext(newprobe, copy_beamline_sptr);
            probecontext.handleAsRing();
            double newhtune, newvtune;
            //newhtune = probecontext.getHorizontalFracTune();
            // newvtune = probecontext.getVerticalFracTune();

            newhtune = probecontext.getHorizontalEigenTune();
            newvtune = probecontext.getVerticalEigenTune();

            dppcount += 1.;
            chromat_H += (newhtune - hcentral_tune) / dpp;
            chromat_V += (newvtune - vcentral_tune) / dpp;
        }

        horizontal_chromaticity = chromat_H / dppcount;
        vertical_chromaticity = chromat_V / dppcount;
        have_chromaticities = true;
    }
}

double
Lattice_simulator::get_horizontal_chromaticity()
{
    get_chromaticities();
    return horizontal_chromaticity;
}

double
Lattice_simulator::get_vertical_chromaticity()
{
    get_chromaticities();
    return vertical_chromaticity;
}

void
write_correctors(Lattice_elements const& horizontal_correctors,
        Lattice_elements const & vertical_correctors,
        Chef_lattice & chef_lattice, std::ofstream & file)
{

    for (Lattice_elements::const_iterator le_it = horizontal_correctors.begin();
            le_it != horizontal_correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin();
                ce_it != chef_elements.end(); ++ce_it) {

            double k2 = 2. * (*ce_it)->Strength() / chef_lattice.get_brho();
            (*le_it)->set_double_attribute("k2", k2);
            file << (*le_it)->get_name() << ":  SEXTUPOLE,  L="
                    << std::setprecision(5)
                    << (*le_it)->get_double_attribute("l") << ",    K2="
                    << std::setprecision(11)
                    << (*le_it)->get_double_attribute("k2") << std::endl;

        }
    }

    for (Lattice_elements::const_iterator le_it = vertical_correctors.begin();
            le_it != vertical_correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin();
                ce_it != chef_elements.end(); ++ce_it) {

            double k2 = 2. * (*ce_it)->Strength() / chef_lattice.get_brho();
            (*le_it)->set_double_attribute("k2", k2);
            file << (*le_it)->get_name() << ":  SEXTUPOLE,  L="
                    << std::setprecision(5)
                    << (*le_it)->get_double_attribute("l") << ",    K2="
                    << std::setprecision(11)
                    << (*le_it)->get_double_attribute("k2") << std::endl;

        }
    }

}

// extract_sextuploe_strengths is a local function
void
extract_sextupole_strengths(Lattice_elements const& correctors,
        Chef_lattice & chef_lattice)
{
    for (Lattice_elements::const_iterator le_it = correctors.begin();
            le_it != correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin();
                ce_it != chef_elements.end(); ++ce_it) {

            double k2 = 2. * (*ce_it)->Strength() / chef_lattice.get_brho();
            (*le_it)->set_double_attribute("k2", k2);

        }
    }
}

void
set_chef_chrom_correctors(Lattice_elements const& correctors,
        Chef_lattice & chef_lattice, BeamlineContext & beamline_context,
        bool horizontal)
{
    for (Lattice_elements::const_iterator le_it = correctors.begin();
            le_it != correctors.end(); ++le_it) {
        Chef_elements chef_elements(chef_lattice.get_chef_elements(*(*le_it)));
        for (Chef_elements::iterator ce_it = chef_elements.begin();
                ce_it != chef_elements.end(); ++ce_it) {

            if ((std::strcmp((*ce_it)->Type(), "sextupole") == 0)
                    && ((!(*le_it)->has_double_attribute("tilt")
                            && !(*le_it)->has_string_attribute("tilt"))
                            || ((*le_it)->has_double_attribute("tilt")
                                    && fabs(
                                            (*le_it)->get_double_attribute(
                                                    "tilt")) < 1.e-6))) {
                if (horizontal) {
                    beamline_context.addHChromCorrector(
                            boost::dynamic_pointer_cast<sextupole >(*ce_it));
                } else {
                    beamline_context.addVChromCorrector(
                            boost::dynamic_pointer_cast<sextupole >(*ce_it));
                }
            } else if ((std::strcmp((*ce_it)->Type(), "thinSextupole") == 0)
                    && ((!(*le_it)->has_double_attribute("tilt")
                            && !(*le_it)->has_string_attribute("tilt"))
                            || ((*le_it)->has_double_attribute("tilt")
                                    && fabs(
                                            (*le_it)->get_double_attribute(
                                                    "tilt")) < 1.e-6))) {

                if (horizontal) {
                    beamline_context.addHChromCorrector(
                            boost::dynamic_pointer_cast<thinSextupole >(
                                    *ce_it));
                } else {
                    beamline_context.addVChromCorrector(
                            boost::dynamic_pointer_cast<thinSextupole >(
                                    *ce_it));
                }
            } else {
                std::string message(
                        "Lattice_simulator::adjust_chromaticities: Lattice_element ");
                message += (*le_it)->get_name();
                message += " of type ";
                message += (*le_it)->get_type();
                message += " cannot be used as a corrector because it has a";
                message += " chef element of type ";
                message += (*ce_it)->Type();
                message += " or it is skewed (i.e. nonzero tilt) ";
                throw std::runtime_error(message.c_str());
            }
        }
    }
}

void
Lattice_simulator::adjust_chromaticities(double horizontal_chromaticity,
        double vertical_chromaticity,
        Lattice_elements const& horizontal_correctors,
        Lattice_elements const& vertical_correctors, double tolerance,
        int max_steps)
{
    if (Jet__environment::getLastEnv() == 0) {
        JetParticle::createStandardEnvironments(map_order);
    }
    BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr());
    BeamlineContext beamline_context(
            reference_particle_to_chef_particle(
                    lattice_sptr->get_reference_particle()), beamline_sptr);
    beamline_context.handleAsRing();
    set_chef_chrom_correctors(horizontal_correctors, *chef_lattice_sptr,
            beamline_context, true);
    set_chef_chrom_correctors(vertical_correctors, *chef_lattice_sptr,
            beamline_context, false);

    double chr_h = get_horizontal_chromaticity();
    double chr_v = get_vertical_chromaticity();
    int rank = Commxx().get_rank();
    if (rank == 0) {
        std::cout << "        Initial chromaticity: horizontal    : " << chr_h
                << std::endl;
        ;
        std::cout << "                       vertical      : " << chr_v
                << std::endl;
        std::cout << "        Desired chromaticity: horizontal      : "
                << horizontal_chromaticity << std::endl;
        std::cout << "                     vertical        : "
                << vertical_chromaticity << std::endl;
    }

    double dh = horizontal_chromaticity - chr_h;
    double dv = vertical_chromaticity - chr_v;
    int count = 0;

    while (((fabs(dh) > tolerance) || (fabs(dv) > tolerance))
            && (count < max_steps)) {

        std::cout << " chr_h=" << chr_h << "  chr_v  =" << chr_v << "  count="
                << count << std::endl;
        std::cout << " dh=" << dh << "  dv=" << dv << "  count=" << count
                << std::endl;
        int status = beamline_context.changeChromaticityBy(dh, dv);

        if (status == BeamlineContext::NO_CHROMATICITY_ADJUSTER) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_chromaticities: no corrector elements found");
        } else if (status != BeamlineContext::OKAY) {
            throw std::runtime_error(
                    "Lattice_simulator::adjust_chromaticities: failed with unknown status");
        }

        extract_sextupole_strengths(horizontal_correctors, *chef_lattice_sptr);
        extract_sextupole_strengths(vertical_correctors, *chef_lattice_sptr);

        have_chromaticities = false;
        chr_h = get_horizontal_chromaticity();
        chr_v = get_vertical_chromaticity();
        dh = horizontal_chromaticity - chr_h;
        dv = vertical_chromaticity - chr_v;
        count++;

    }

    if (rank == 0) {
        std::ofstream file;
        file.open("sextupole_correctors.txt");
        write_correctors(horizontal_correctors, vertical_correctors,
                *chef_lattice_sptr, file);
        file.close();
    }

    have_chromaticities = false;
    if (rank == 0) {
        std::cout << " Chromaticity adjusted in " << count << " steps"
                << std::endl;
        if (count == max_steps) std::cout
                << " Convergence not reached, increase the maximum number of steps"
                << std::endl;
        std::cout << "  final    chromaticity:  horizontal    : " << chr_h
                << "     vertical     : " << chr_v << std::endl;
    }

}

template<class Archive>
    void
    Lattice_simulator::save(Archive & ar, const unsigned int version) const
    {
        ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
        ar & BOOST_SERIALIZATION_NVP(slices);
        ar & BOOST_SERIALIZATION_NVP(have_slices);
        ar & BOOST_SERIALIZATION_NVP(chef_lattice_sptr);
        ar & BOOST_SERIALIZATION_NVP(extractor_map_sptr);
        ar & BOOST_SERIALIZATION_NVP(aperture_extractor_map_sptr);
        ar & BOOST_SERIALIZATION_NVP(have_beamline_context);
        ar & BOOST_SERIALIZATION_NVP(map_order);
        ar & BOOST_SERIALIZATION_NVP(bucket_length);
        ar & BOOST_SERIALIZATION_NVP(have_element_lattice_functions);
        ar & BOOST_SERIALIZATION_NVP(have_slice_lattice_functions);
        ar & BOOST_SERIALIZATION_NVP(horizontal_tune);
        ar & BOOST_SERIALIZATION_NVP(vertical_tune);
        ar & BOOST_SERIALIZATION_NVP(have_tunes);
        ar & BOOST_SERIALIZATION_NVP(horizontal_chromaticity);
        ar & BOOST_SERIALIZATION_NVP(vertical_chromaticity);
        ar & BOOST_SERIALIZATION_NVP(have_chromaticities);
        ar & BOOST_SERIALIZATION_NVP(lattice_functions_element_map);
        ar & BOOST_SERIALIZATION_NVP(lattice_functions_slice_map);
        ar & BOOST_SERIALIZATION_NVP(linear_one_turn_map);
    }
template<class Archive>
    void
    Lattice_simulator::load(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(lattice_sptr);
        ar & BOOST_SERIALIZATION_NVP(slices);
        ar & BOOST_SERIALIZATION_NVP(have_slices);
        ar & BOOST_SERIALIZATION_NVP(chef_lattice_sptr);
        ar & BOOST_SERIALIZATION_NVP(extractor_map_sptr);
        ar & BOOST_SERIALIZATION_NVP(aperture_extractor_map_sptr);
        ar & BOOST_SERIALIZATION_NVP(have_beamline_context);
        ar & BOOST_SERIALIZATION_NVP(map_order);
        ar & BOOST_SERIALIZATION_NVP(bucket_length);
        ar & BOOST_SERIALIZATION_NVP(have_element_lattice_functions);
        ar & BOOST_SERIALIZATION_NVP(have_slice_lattice_functions);
        ar & BOOST_SERIALIZATION_NVP(horizontal_tune);
        ar & BOOST_SERIALIZATION_NVP(vertical_tune);
        ar & BOOST_SERIALIZATION_NVP(have_tunes);
        ar & BOOST_SERIALIZATION_NVP(horizontal_chromaticity);
        ar & BOOST_SERIALIZATION_NVP(vertical_chromaticity);
        ar & BOOST_SERIALIZATION_NVP(have_chromaticities);
        ar & BOOST_SERIALIZATION_NVP(lattice_functions_element_map);
        ar & BOOST_SERIALIZATION_NVP(lattice_functions_slice_map);
        ar & BOOST_SERIALIZATION_NVP(linear_one_turn_map);
        if (have_beamline_context) {
            calculate_beamline_context();
        }
        normal_form_sage_sptr.reset();
    }

template
void
Lattice_simulator::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Lattice_simulator::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Lattice_simulator::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Lattice_simulator::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Lattice_simulator::~Lattice_simulator()
{
}

