#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/propagator.h"
#include "components/simulation/lattice_simulator.h"
#include "components/foundation/physical_constants.h"
#include "components/lattice/chef_utils.h"
#include "components/bunch/bunch.h"
#include "lattice_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 7);

    Propagator propagator(stepper);
}

BOOST_FIXTURE_TEST_CASE(propagate, Lattice_fixture)
{
    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 4);
    Propagator propagator(stepper);

    int num_turns = 4;
    propagator.propagate(b.bunch, num_turns, true, true);
}

