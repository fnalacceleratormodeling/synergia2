#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/lattice_simulator.h"
#include "lattice_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
}

