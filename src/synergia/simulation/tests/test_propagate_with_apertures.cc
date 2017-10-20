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

struct Big_fixture
{
    Big_fixture() : l(), b()
    {
        BOOST_TEST_MESSAGE("setup big fixture");
    }
    ~Big_fixture()
    {
        BOOST_TEST_MESSAGE("teardown big fixture");
    }
    Lattice_fixture l;
    Bunch_fixture b;
}
;
        
int
particle_mask(MArray2d_ref m, int np)
{
    int mask = 0;
    for (int i=0; i<np; ++i) {
        int partnum = m[i][Bunch::id];
        mask |= (1 << partnum);
    }
    return mask;
}

// independent stepper 1 step
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_aperture_1, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper(l.lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper 2 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_aperture_2, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper(l.lattice_sptr, 1, 2));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper 3 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_aperture_3, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper(l.lattice_sptr, 1, 3));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper 7 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_aperture_7, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper(l.lattice_sptr, 1, 7));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper elements 1 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_elements_aperture_1, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(l.lattice_sptr, 1, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper elements 2 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_elements_aperture_2, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(l.lattice_sptr, 1, 2));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper elements 3 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_elements_aperture_3, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(l.lattice_sptr, 1, 3));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// independent stepper elements 7 steps
BOOST_FIXTURE_TEST_CASE(propagate_independent_stepper_elements_aperture_7, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Stepper_sptr stepper_sptr(new Independent_stepper_elements(l.lattice_sptr, 1, 7));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper 1 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_aperture_1, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper(l.lattice_sptr, 1, dummy, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper 2 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_aperture_2, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper(l.lattice_sptr, 1, dummy, 2));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper 3 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_aperture_3, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper(l.lattice_sptr, 1, dummy, 3));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper 7 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_aperture_7, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper(l.lattice_sptr, 1, dummy, 7));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper elements 1 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_elements_aperture_1, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper_elements(l.lattice_sptr, 1, dummy, 1));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper elements 2 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_elements_aperture_2, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper_elements(l.lattice_sptr, 1, dummy, 2));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper elements 3 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_elements_aperture_3, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper_elements(l.lattice_sptr, 1, dummy, 3));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}

// split operator stepper elements 7 steps
BOOST_FIXTURE_TEST_CASE(propagate_splitoperator_stepper_elements_aperture_7, Big_fixture)
{
    Bunch_simulator bunch_simulator(b.bunch_sptr);
    Collective_operator_sptr dummy(new Dummy_collective_operator("foo"));
    Stepper_sptr stepper_sptr(new Split_operator_stepper_elements(l.lattice_sptr, 1, dummy, 7));
    Propagator propagator(stepper_sptr);
    propagator.propagate(bunch_simulator, 1, 1, verbosity);
    // check particles 0 and 2 survive
    BOOST_CHECK_EQUAL(particle_mask(b.bunch_sptr->get_local_particles(), b.bunch_sptr->get_local_num()), 5);
}
