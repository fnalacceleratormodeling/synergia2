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

double tolerance = 1.0e-13;

const double quad_length = 0.2;
const double drift_length = 0.8;
const double k1 = 1.0/(quad_length * 0.7);

const std::string name("foo_lattice");

// propagate a test particle through a lattice
MArray1d
propagate_lattice(Lattice_sptr lattice_sptr, MArray1d_ref co)
{
    Commxx_sptr commxx(new Commxx());
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(), commxx->get_size(), 1.0e10, commxx));
    for (int i=0; i<6; ++i) {
        bunch_sptr->get_local_particles()[0][i] = co[i];
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    Independent_stepper_sptr stepper_sptr(new Independent_stepper(lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, 0);
    
    MArray1d results(boost::extents[6]);
    for (int i=0; i<6; ++i) {
        results[i] = bunch_sptr->get_local_particles()[0][i];
    }
    return results;
}

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
    MArray1d result(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
                     closed_orbit[i], co_tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(ring_on_momentum)
{
    Lsexpr foborodobo32lsx(read_lsexpr_file("lattices/foborodobo32.lsx"));
    Lattice_sptr lattice_sptr(new Lattice(foborodobo32lsx));

    MArray1d closed_orbit(boost::extents[6]);
    MArray1d result(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
                     closed_orbit[i], co_tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(ring_off_momentum)
{
    Lsexpr foborodobo32lsx(read_lsexpr_file("lattices/foborodobo32.lsx"));
    Lattice_sptr lattice_sptr(new Lattice(foborodobo32lsx));

    MArray1d closed_orbit(boost::extents[6]);
    MArray1d result(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
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
    MArray1d result(boost::extents[6]);

    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
                     closed_orbit[i], co_tolerance) );
    }
}
#endif

BOOST_AUTO_TEST_CASE(closed_orbit_with_kicker)
{
    // read the lattice
    Lattice_sptr lattice_sptr(MadX_reader().get_lattice_sptr("model", "lattices/foborodobo128.madx"));
    // turn on one kicker
    Lattice_elements elements(lattice_sptr->get_elements());
    bool found_kicker = false;
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        // first kicker is hc1
        if ((*lit)->get_name() == "hc1") {
            (*lit)->set_double_attribute("kick", 0.02);
            found_kicker = true;
            break;
        }
    }
    if (! found_kicker) {
        throw std::runtime_error("did not find kicker in foborodobo128 lattice");
    }

    MArray1d closed_orbit(boost::extents[6]);
    MArray1d result(boost::extents[6]);

    // 1.0e-13 works, 1.0e-14 doesn't converge within 30 iterations
    closed_orbit = calculate_closed_orbit(lattice_sptr, 0.0, 1.0e-13);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
                     closed_orbit[i], co_tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(closed_orbit_off_momentum_with_kicker)
{
    // read the lattice
    Lattice_sptr lattice_sptr(MadX_reader().get_lattice_sptr("model", "lattices/foborodobo128.madx"));
    // turn on one kicker
    Lattice_elements elements(lattice_sptr->get_elements());
    bool found_kicker = false;
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        // first kicker is hc1
        if ((*lit)->get_name() == "hc1") {
            (*lit)->set_double_attribute("kick", 0.02);
            found_kicker = true;
            break;
        }
    }
    if (! found_kicker) {
        throw std::runtime_error("did not find kicker in foborodobo128 lattice");
    }

    MArray1d closed_orbit(boost::extents[6]);
    MArray1d result(boost::extents[6]);

    // 1.0e-13 works, 1.0e-14 doesn't converge within 30 iterations
    closed_orbit = calculate_closed_orbit(lattice_sptr, 1.0e-2, 1.0e-13);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
                     closed_orbit[i], co_tolerance) );
    }
}

BOOST_AUTO_TEST_CASE(closed_orbit_off_momentum_with_kicker_and_skew_element)
{
    // read the lattice
    Lattice_sptr lattice_sptr(MadX_reader().get_lattice_sptr("model", "lattices/foborodobo128.madx"));
    // turn on one kicker
    Lattice_elements elements(lattice_sptr->get_elements());
    bool found_kicker = false;
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        // first kicker is hc1
        if ((*lit)->get_name() == "hc1") {
            (*lit)->set_double_attribute("kick", 0.02);
            found_kicker = true;
            break;
        }
    }
    if (! found_kicker) {
        throw std::runtime_error("did not find kicker in foborodobo128 lattice");
    }
    for (Lattice_elements::const_iterator lit=elements.begin(); lit!=elements.end(); ++lit) {
        // add skew component to the first focussing quadrupole
        if (((*lit)->get_type() == "quadrupole") &&
             ((*lit)->get_double_attribute("k1") > 0.0)) {
            double k1 = (*lit)->get_double_attribute("k1");
            double k1s = k1 * std::sin(mconstants::pi/48.0);
            k1 *= std::cos(mconstants::pi/48.0);
            break;
        }
    }

    MArray1d closed_orbit(boost::extents[6]);
    MArray1d result(boost::extents[6]);

    // 1.0e-13 works, 1.0e-14 doesn't converge within 30 iterations
    closed_orbit = calculate_closed_orbit(lattice_sptr, 1.0e-2, 1.0e-13);
    // Now let's see if it actually works

    result = propagate_lattice(lattice_sptr, closed_orbit);

    const double co_tolerance = 100.0 * tolerance * lattice_sptr->get_length();

    for (int i=0; i<4; ++i) {
        BOOST_CHECK( floating_point_equal(result[i],
                     closed_orbit[i], co_tolerance) );
    }
}

