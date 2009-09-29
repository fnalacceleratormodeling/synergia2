#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include "parallel_utils.h"

struct MPI_fixture
{
    MPI_fixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        MPI_Init(&argc, &argv);
    }
    ~MPI_fixture()
    {
        MPI_Finalize();
    }
};

BOOST_GLOBAL_FIXTURE(MPI_fixture);


BOOST_AUTO_TEST_CASE(test_mpi_get_rank)
{
    int rank = mpi_get_rank(MPI_COMM_WORLD);
    BOOST_CHECK_EQUAL(rank,0);
}
