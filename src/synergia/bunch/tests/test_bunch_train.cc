#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/train.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
}

BOOST_AUTO_TEST_CASE(get_num_bunches)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    BOOST_CHECK_EQUAL(bunch_train.get_num_bunches(), num_bunches);
}

BOOST_AUTO_TEST_CASE(get_bunch_separation)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    BOOST_CHECK_EQUAL(bunch_train.get_bunch_separation(), bunch_separation);
}

BOOST_AUTO_TEST_CASE(get_master_comm)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    BOOST_CHECK_EQUAL(bunch_train.get_master_comm().get(), MPI_COMM_WORLD);
}

BOOST_AUTO_TEST_CASE(get_comm)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        bunch_train.get_comm(bunch);
        // don't know what to check here...
    }
}

BOOST_AUTO_TEST_CASE(get_comm_bad_index)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    bool caught = false;
    try {
        bunch_train.get_comm(num_bunches + 1);
    } catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(is_on_this_rank)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        BOOST_CHECK(bunch_train.is_on_this_rank(bunch));
    }
}

BOOST_AUTO_TEST_CASE(get_bunch_sptr)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    for (int bunch = 0; bunch < num_bunches; ++bunch) {
        bunch_train.get_bunch_sptr(bunch);
        // don't know what to check here...
    }
}

BOOST_AUTO_TEST_CASE(get_bunch_sptr_bad_index)
{
    const int num_bunches = 3;
    const double bunch_separation = 1.7;
    Commxx_sptr commxx_sptr(new Commxx);
    Bunch_train bunch_train(num_bunches, bunch_separation,
            commxx_sptr);
    bool caught = false;
    try {
        bunch_train.get_bunch_sptr(num_bunches + 1);
    } catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}


