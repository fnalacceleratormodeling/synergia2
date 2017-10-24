#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture);

const int verbosity = 1;
const double tolerance = 1.0e-12;

// fixture includes a lattice of marker, drift, drift with apertures on the first
// two elements.  A bunch contains three particles, two will make it through and one
// will be stopped by the aperture.

const double total_energy = 4.0;
const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.0e11;
const int local_num = 3;

const double driftlen = 5.0;
const double aperture_radius = 0.01;

const double eps = 1.0e-10;

struct Reference_particle_fixture
{
    Reference_particle_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge, four_momentum)
    {
        BOOST_TEST_MESSAGE("setup reference particle fixture");
    }
    ~Reference_particle_fixture()
    {
        BOOST_TEST_MESSAGE("teardown reference particle fixture");
    }
    Four_momentum four_momentum;
    Reference_particle reference_particle;
}
;

struct Lattice_fixture
{
        Lattice_fixture() :
            rpf(), lattice_sptr(new Lattice())
        {
                BOOST_TEST_MESSAGE("setup lattice fixture");

                // aperture first two elements
        Lattice_element mrk("marker", "mrk");
        mrk.set_string_attribute("aperture_type", "circular");
        mrk.set_double_attribute("circular_aperture_radius", aperture_radius);

        Lattice_element drft1("drift", "bar");
        drft1.set_string_attribute("aperture_type", "circular");
        drft1.set_double_attribute("circular_aperture_radius", aperture_radius);
        drft1.set_double_attribute("l", driftlen);

                // no aperture third element
        Lattice_element drft2("drift", "baz");
        drft2.set_double_attribute("l", driftlen);

                lattice_sptr->append(mrk);
                lattice_sptr->append(drft1);
                lattice_sptr->append(drft2);

        lattice_sptr->set_reference_particle(rpf.reference_particle);
    }

    ~Lattice_fixture()
    {
        BOOST_TEST_MESSAGE("teardown lattice fixture");
    }
    Reference_particle_fixture rpf;
    Lattice_sptr lattice_sptr;
}
;

struct Bunch_fixture
{
    Bunch_fixture() : rpf(), comm_sptr(new Commxx), bunch_sptr()
    {
        BOOST_TEST_MESSAGE("setup bunch fixture");
        bunch_sptr = Bunch_sptr(new Bunch(rpf.reference_particle, local_num, real_num, comm_sptr));
        MArray2d_ref local_particles(bunch_sptr->get_local_particles());
        // three particles: the first goes through the center of the lattice
        // the second has enough just transverse momentum to exit the aperture of the first drift
        // the third has just barely misses the aperture cutoff
        local_particles[1][1] = aperture_radius / std::sqrt(aperture_radius*aperture_radius + driftlen*driftlen) + eps;
        local_particles[2][3] = aperture_radius / std::sqrt(aperture_radius*aperture_radius + driftlen*driftlen) - eps;
        // the particle numbering gets messed up when bunch is called over and over again
        for (int i=0; i<3; ++i) {
            local_particles[i][Bunch::id] = i;
        }
    }
    ~Bunch_fixture()
    {
        BOOST_TEST_MESSAGE("teardown bunch fixture");
    }
    Reference_particle_fixture rpf;
    Commxx_sptr comm_sptr;
    Bunch_sptr bunch_sptr;
}
;


// test that setting extractor_type to a linear chef_map works for a nonzero length sextupole
// 1 step
BOOST_AUTO_TEST_CASE(test_sextupole_chef_map_extractor1)
{
    const double sexlength = 1.0;
    const double sexstrength = 5.0;

    Reference_particle_fixture rpf;

    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element mrk("marker", "mrk1");
    Lattice_element drft1("drift", "d1");
    drft1.set_double_attribute("l", driftlen);

    Lattice_element sex("sextupole", "s");
    sex.set_double_attribute("l", sexlength);
    sex.set_double_attribute("k2", sexstrength);
    sex.set_string_attribute("extractor_type", "chef_map");

    lattice_sptr->append(mrk);
    lattice_sptr->append(drft1);
    lattice_sptr->append(mrk);
    lattice_sptr->append(sex);
    lattice_sptr->append(drft1);

    lattice_sptr->set_reference_particle(rpf.reference_particle);

    Commxx_sptr comm(new Commxx);

    Bunch_sptr bunch_sptr(new Bunch(rpf.reference_particle, 1, real_num, comm));
    Bunch_simulator bunch_simulator(bunch_sptr);
    MArray2d_ref local_particles(bunch_sptr->get_local_particles());
    local_particles[0][Bunch::x] = 0.013;
    local_particles[0][Bunch::y] = 0.012;

    Stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));

    Propagator propagator(stepper_sptr);

    propagator.propagate(bunch_simulator, 1, 1, 1);

#if 0
    for (int i=0; i<6; ++i) {
        std::cout << i << ": " << local_particles[0][i] << std::endl;
    }
#endif

    // If chef_maps is working, propagating through the sextupole with linear maps
    // should not affect transverse momentum.

    BOOST_CHECK_EQUAL(local_particles[0][Bunch::x], 0.013);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::xp], 0.0);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::y], 0.012);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::yp], 0.0);

}

// test that setting extractor_type to a linear chef_map works for a nonzero length sextupole
// 1 step/element
BOOST_AUTO_TEST_CASE(test_sextupole_chef_map_extractor1e)
{
    const double sexlength = 1.0;
    const double sexstrength = 5.0;

    Reference_particle_fixture rpf;

    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element mrk("marker", "mrk1");
    Lattice_element drft1("drift", "d1");
    drft1.set_double_attribute("l", driftlen);

    Lattice_element sex("sextupole", "s");
    sex.set_double_attribute("l", sexlength);
    sex.set_double_attribute("k2", sexstrength);
    sex.set_string_attribute("extractor_type", "chef_map");

    lattice_sptr->append(mrk);
    lattice_sptr->append(drft1);
    lattice_sptr->append(mrk);
    lattice_sptr->append(sex);
    lattice_sptr->append(drft1);

    lattice_sptr->set_reference_particle(rpf.reference_particle);

    Commxx_sptr comm(new Commxx);

    Bunch_sptr bunch_sptr(new Bunch(rpf.reference_particle, 1, real_num, comm));
    Bunch_simulator bunch_simulator(bunch_sptr);
    MArray2d_ref local_particles(bunch_sptr->get_local_particles());
    local_particles[0][Bunch::x] = 0.013;
    local_particles[0][Bunch::y] = 0.012;

    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 1));

    Propagator propagator(stepper_sptr);

    propagator.propagate(bunch_simulator, 1, 1, 1);

#if 0
    for (int i=0; i<6; ++i) {
        std::cout << i << ": " << local_particles[0][i] << std::endl;
    }
#endif

    // If chef_maps is working, propagating through the sextupole with linear maps
    // should not affect transverse momentum.

    BOOST_CHECK_EQUAL(local_particles[0][Bunch::x], 0.013);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::xp], 0.0);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::y], 0.012);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::yp], 0.0);

}

// test that setting extractor_type to a linear chef_map works for a nonzero length sextupole
// 8 step
BOOST_AUTO_TEST_CASE(test_sextupole_chef_map_extractor8)
{
    const double sexlength = 1.0;
    const double sexstrength = 5.0;

    Reference_particle_fixture rpf;

    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element mrk("marker", "mrk1");
    Lattice_element drft1("drift", "d1");
    drft1.set_double_attribute("l", driftlen);

    Lattice_element sex("sextupole", "s");
    sex.set_double_attribute("l", sexlength);
    sex.set_double_attribute("k2", sexstrength);
    sex.set_string_attribute("extractor_type", "chef_map");

    lattice_sptr->append(mrk);
    lattice_sptr->append(drft1);
    lattice_sptr->append(mrk);
    lattice_sptr->append(sex);
    lattice_sptr->append(drft1);

    lattice_sptr->set_reference_particle(rpf.reference_particle);

    Commxx_sptr comm(new Commxx);

    Bunch_sptr bunch_sptr(new Bunch(rpf.reference_particle, 1, real_num, comm));
    Bunch_simulator bunch_simulator(bunch_sptr);
    MArray2d_ref local_particles(bunch_sptr->get_local_particles());
    local_particles[0][Bunch::x] = 0.013;
    local_particles[0][Bunch::y] = 0.012;

    Stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 8));

    Propagator propagator(stepper_sptr);

    propagator.propagate(bunch_simulator, 1, 8, 1);

#if 0
    for (int i=0; i<6; ++i) {
        std::cout << i << ": " << local_particles[0][i] << std::endl;
    }
#endif
    // If chef_maps is working, propagating through the sextupole with linear maps
    // should not affect transverse momentum.

    BOOST_CHECK_EQUAL(local_particles[0][Bunch::x], 0.013);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::xp], 0.0);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::y], 0.012);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::yp], 0.0);

}

// test that setting extractor_type to a linear chef_map works for a nonzero length sextupole
// 3 step/ement
BOOST_AUTO_TEST_CASE(test_sextupole_chef_map_extractor3e)
{
    const double sexlength = 1.0;
    const double sexstrength = 5.0;

    Reference_particle_fixture rpf;

    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element mrk("marker", "mrk1");
    Lattice_element drft1("drift", "d1");
    drft1.set_double_attribute("l", driftlen);

    Lattice_element sex("sextupole", "s");
    sex.set_double_attribute("l", sexlength);
    sex.set_double_attribute("k2", sexstrength);
    sex.set_string_attribute("extractor_type", "chef_map");

    lattice_sptr->append(mrk);
    lattice_sptr->append(drft1);
    lattice_sptr->append(mrk);
    lattice_sptr->append(sex);
    lattice_sptr->append(drft1);

    lattice_sptr->set_reference_particle(rpf.reference_particle);

    Commxx_sptr comm(new Commxx);

    Bunch_sptr bunch_sptr(new Bunch(rpf.reference_particle, 1, real_num, comm));
    Bunch_simulator bunch_simulator(bunch_sptr);
    MArray2d_ref local_particles(bunch_sptr->get_local_particles());
    local_particles[0][Bunch::x] = 0.013;
    local_particles[0][Bunch::y] = 0.012;

    Stepper_sptr stepper_sptr(new Independent_stepper_elements(lattice_sptr, 1, 3));

    Propagator propagator(stepper_sptr);

    propagator.propagate(bunch_simulator, 1, 3, 1);

#if 0
    for (int i=0; i<6; ++i) {
        std::cout << i << ": " << local_particles[0][i] << std::endl;
    }
#endif
    // If chef_maps is working, propagating through the sextupole with linear maps
    // should not affect transverse momentum.

    BOOST_CHECK_EQUAL(local_particles[0][Bunch::x], 0.013);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::xp], 0.0);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::y], 0.012);
    BOOST_CHECK_EQUAL(local_particles[0][Bunch::yp], 0.0);

}

BOOST_AUTO_TEST_CASE(test_rfcavity_chef_map_extractor1)
{
    const double rfvoltage = 0.05;

    Reference_particle_fixture rpf;

    Lattice_sptr lattice_sptr(new Lattice("foo"));

    Lattice_element mrk("marker", "mrk1");
    Lattice_element drft1("drift", "d1");
    drft1.set_double_attribute("l", driftlen);

    Lattice_element cav("rfcavity", "r");
    cav.set_double_attribute("l", 0.0);
    cav.set_double_attribute("volt", rfvoltage);
    cav.set_double_attribute("harmon", 1.0);
    cav.set_string_attribute("extractor_type", "chef_map");

    lattice_sptr->append(mrk);
    lattice_sptr->append(drft1);
    lattice_sptr->append(mrk);
    lattice_sptr->append(cav);

    lattice_sptr->set_reference_particle(rpf.reference_particle);

    double beta = rpf.reference_particle.get_beta();
    double pmom = rpf.reference_particle.get_momentum();

    double lattice_length = lattice_sptr->get_length();
    double omega = 2.0 * mconstants::pi * beta * pconstants::c/lattice_length;

    Commxx_sptr comm(new Commxx);

    Bunch_sptr bunch_sptr(new Bunch(rpf.reference_particle, 1, real_num, comm));
    Bunch_simulator bunch_simulator(bunch_sptr);
    MArray2d_ref local_particles(bunch_sptr->get_local_particles());
    // this is exactly one cavity wavelength.  The full cavity propagator would give 0 kick
    local_particles[0][Bunch::cdt] = lattice_length/beta;

    Stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));

    Propagator propagator(stepper_sptr);

    propagator.propagate(bunch_simulator, 1, 1, 1);

#if 0
    for (int i=0; i<6; ++i) {
        std::cout << i << ": " << local_particles[0][i] << std::endl;
    }
#endif

    // the linearized energy change = V * omega/c *  c dt
    double dE = 1.0e-3*rfvoltage * (omega/pconstants::c) * (lattice_length/beta);
    double dp = dE/beta;
    double dpop =dp/pmom;
    BOOST_CHECK_CLOSE(local_particles[0][Bunch::dpop],  dpop, 1.0e-8);

}
