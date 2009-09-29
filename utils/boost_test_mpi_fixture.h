#ifndef BOOST_TEST_MPI_FIXTURE_H_
#define BOOST_TEST_MPI_FIXTURE_H_

#include "mpi.h"

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

#endif /* BOOST_TEST_MPI_FIXTURE_H_ */
