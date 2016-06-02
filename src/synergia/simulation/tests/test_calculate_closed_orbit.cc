#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/floating_point.h"

#include "synergia/simulation/calculate_closed_orbit.h"

double tolerance = 1.0e-15;

const double quad_length = 0.2;
const double drift_length = 0.8;
const double k1 = 1.0/(quad_length * 0.7);

const std::string name("foo_lattice");

BOOST_AUTO_TEST_CASE(trivial_closed_orbit)
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    f.set_double_attribute("k1", k1);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", -k1);

    Lattice_sptr lattice_sptr(new Lattice(name));
    lattice_sptr->append(f);
    lattice_sptr->append(o);
    lattice_sptr->append(d);
    lattice_sptr->append(o);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice_sptr->set_reference_particle(reference_particle);

    MArray1d closed_orbit(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    Commxx_sptr commxx(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(), commxx->get_size(), 1.0e10, commxx));
    for (int i=0; i<6; ++i) {
        bunch_sptr->get_local_particles()[0][i] = closed_orbit[i];
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, 0);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(bunch_sptr->get_local_particles()[0][i],
                     closed_orbit[i], co_tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(ring_on_momentum)
{
    Lsexpr foborodobo32lsx(read_lsexpr_file("lattices/foborodobo32.lsx"));
    Lattice_sptr lattice_sptr(new Lattice(foborodobo32lsx));

    MArray1d closed_orbit(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    Commxx_sptr commxx(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(), commxx->get_size(), 1.0e10, commxx));
    for (int i=0; i<6; ++i) {
        bunch_sptr->get_local_particles()[0][i] = closed_orbit[i];
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, 0);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(bunch_sptr->get_local_particles()[0][i],
                     closed_orbit[i], co_tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(ring_off_momentum)
{
    Lsexpr foborodobo32lsx(read_lsexpr_file("lattices/foborodobo32.lsx"));
    Lattice_sptr lattice_sptr(new Lattice(foborodobo32lsx));

    MArray1d closed_orbit(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.001);
    // Now let's see if it actually works

    Commxx_sptr commxx(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(), commxx->get_size(), 1.0e10, commxx));
    for (int i=0; i<6; ++i) {
        bunch_sptr->get_local_particles()[0][i] = closed_orbit[i];
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, 0);
    
    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(bunch_sptr->get_local_particles()[0][i],
                     closed_orbit[i], co_tolerance) );
    }
}

#if 0
BOOST_AUTO_TEST_CASE(unstable_lattice_has_no_closed_orbit)
{
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    // not focussed enough to be stable
    f.set_double_attribute("k1", k1/10.0);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", -k1/10.0);

    Lattice_sptr lattice_sptr(new Lattice(name));
    lattice_sptr->append(f);
    lattice_sptr->append(o);
    lattice_sptr->append(d);
    lattice_sptr->append(o);
    Reference_particle reference_particle(charge, mass, total_energy);
    lattice_sptr->set_reference_particle(reference_particle);

    MArray1d closed_orbit(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    Commxx_sptr commxx(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(), 1, 1.0e10, commxx));
    for (int i=0; i<6; ++i) {
        bunch_sptr->get_local_particles()[0][i] = closed_orbit[i];
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, 0);

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(bunch_sptr->get_local_particles()[0][i],
                     closed_orbit[i], tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(closed_orbit_with_kicker)
{
#if 0
    // read the lattice
    Lattice_sptr lattice_sptr(MadX_reader().get_lattice_sptr("model", "lattices/foborodobo128.madx"));
    // turn on one kicker
    Lattice_elements elements(lattice_sptr->get_elements());
    bool found_kicker = false;
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        // first kicker is hc1
        if ((*lit)->get_name() == "hc1") {
            //(*lit)->set_double_attribute("kick", 0.2);
            found_kicker = true;
            break;
        }
    }
    if (! found_kicker) {
        throw std::runtime_error("did not find kicker in foborodobo128 lattice");
    }
#endif
    Lattice_sptr orig_lattice_sptr(MadX_reader().get_lattice_sptr("fodo", "../../../../examples/envelope/fodo.seq"));
    Lattice_sptr lattice_sptr(new Lattice(*orig_lattice_sptr));
    Independent_stepper_sptr test_stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator test_propagator(test_stepper_sptr);

    // Now get closed orbit
    MArray1d closed_orbit(boost::extents[6]);
#if 0
    Lattice_sptr lattice_sptr2;
    lattice_sptr2.reset(new Lattice("foo"));
    lattice_sptr2->set_reference_particle(lattice_sptr->get_reference_particle());
    lattice_sptr = lattice_sptr2;
#endif
    std::cout << "Read lattice: " << lattice_sptr->get_name() << " length: " << lattice_sptr->get_length() << ", number elements: " << lattice_sptr->get_elements().size() << std::endl;
    std::cout << "beam energy: " << lattice_sptr->get_reference_particle().get_total_energy() << std::endl;
    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    Commxx_sptr commxx(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(), 1, 1.0e10, commxx));
    for (int i=0; i<6; ++i) {
        bunch_sptr->get_local_particles()[0][i] = closed_orbit[i];
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, 0);

    for (int i=0; i<6; ++i) {
        BOOST_CHECK( floating_point_equal(bunch_sptr->get_local_particles()[0][i],
                     closed_orbit[i], tolerance) );
    }
}
#endif
