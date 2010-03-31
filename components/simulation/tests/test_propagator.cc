#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/propagator.h"
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

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);
    lattice.set_reference_particle(reference_particle);

    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));

    Split_operator_stepper stepper(lattice, 7, space_charge);

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

    Lattice lattice(name);
    lattice.append(f);
    lattice.append(o);
    lattice.append(d);
    lattice.append(o);
    lattice.set_reference_particle(reference_particle);

    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));

    Split_operator_stepper stepper(lattice, 4, space_charge);
    Propagator propagator(stepper);

    Commxx comm(MPI_COMM_WORLD);
    Bunch bunch(reference_particle, constants::proton_charge, total_num, real_num, comm);

    int num_turns = 4;
    propagator.propagate(bunch, num_turns, true, true);
}

