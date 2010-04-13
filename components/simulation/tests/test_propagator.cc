#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/propagator.h"
#include "components/simulation/lattice_simulator.h"
#include "components/foundation/physical_constants.h"
#include "components/lattice/chef_utils.h"
#include "components/bunch/bunch.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const std::string name("foo");
const double mass = constants::mp;
const double real_num = 1.0e11;
const int total_num = 1000;
const double total_energy = 125.0;
const double tolerance = 1.0e-12;
const double quad_length = 0.2;
const double drift_length = 3.0;
const double bend_length = 4.0;

BOOST_AUTO_TEST_CASE(construct)
{

    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice_sptr lattice_sptr(new Lattice(name));
    lattice_sptr->append(f);
    lattice_sptr->append(o);
    lattice_sptr->append(d);
    lattice_sptr->append(o);
    lattice_sptr->set_reference_particle(reference_particle);

    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 7);

    Propagator propagator(stepper);
}

BOOST_AUTO_TEST_CASE(propagate)
{

    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);
    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);

    Lattice_sptr lattice_sptr(new Lattice(name));
    lattice_sptr->append(f);
    lattice_sptr->append(o);
    lattice_sptr->append(d);
    lattice_sptr->append(o);
    lattice_sptr->set_reference_particle(reference_particle);

    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 4);
    Propagator propagator(stepper);

    Commxx comm(MPI_COMM_WORLD);
    Bunch bunch(reference_particle, constants::proton_charge, total_num,
            real_num, comm);

    int num_turns = 4;
    propagator.propagate(bunch, num_turns, true, true);
}

