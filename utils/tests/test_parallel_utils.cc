#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/boost_test_mpi_fixture.h"
#include "parallel_utils.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture);

BOOST_AUTO_TEST_CASE(test_mpi_get_rank)
{
    int rank = mpi_get_rank(MPI_COMM_WORLD);
    BOOST_CHECK_EQUAL(rank,0);
}
