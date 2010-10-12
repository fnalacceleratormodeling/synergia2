#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/collective/space_charge_3d_open_hockney.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct1)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx comm(MPI_COMM_WORLD);

    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
}
