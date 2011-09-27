#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/bunch_train.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct1)
{
    const int num_bunches = 1;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
}

BOOST_AUTO_TEST_CASE(construct2)
{
    const int num_bunches = 2;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
}

BOOST_AUTO_TEST_CASE(construct3)
{
    const int num_bunches = 16;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
}

BOOST_AUTO_TEST_CASE(get_master_comm)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
    BOOST_CHECK_EQUAL(bunch_train.get_master_comm().get(), MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(get_comm)
{
    const int num_bunches = 16;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        bunch_train.get_comm(bunch);
        // don't know what to check here...
    }
}

BOOST_AUTO_TEST_CASE(is_on_this_rank)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(num_bunches, bunch_separation,
            Commxx(MPI_COMM_WORLD));
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        bool on_this_rank = bunch_train.is_on_this_rank(bunch);
        // could use a better test here...
    }
}
