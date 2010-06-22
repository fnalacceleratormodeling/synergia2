#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/distributed_fft3d.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

// n.b. We use 0,1,2 here instead of x,y,z because
// we may use z,y,x ordering of arrays.
const int shape0 = 3;
const int shape1 = 4;
const int shape2 = 5;
const bool periodic = false;

BOOST_AUTO_TEST_CASE(construct)
{
    std::vector<int > shape(3);
    shape[0] = shape0;
    shape[1] = shape1;
    shape[2] = shape2;
    Distributed_fft3d distrubuted_fft3d(shape, periodic, Commxx(MPI_COMM_WORLD));
}

