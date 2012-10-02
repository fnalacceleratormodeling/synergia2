#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/bunch_train.h"
#include "bunches_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_FIXTURE_TEST_CASE(construct1, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);
}

BOOST_FIXTURE_TEST_CASE(construct2, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);
}

BOOST_FIXTURE_TEST_CASE(construct3, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    const double bunch_separation3 = 99.9;
    separations.push_back(bunch_separation3);
    bool caught_error = false;
    try {
        Bunch_train bunch_train(bunches, separations);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(get_num_bunches, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);

    BOOST_CHECK_EQUAL(bunch_train.get_size(), num_bunches);
}

BOOST_FIXTURE_TEST_CASE(get_bunches, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);

    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(0), bunches.at(0));
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(1), bunches.at(1));
    BOOST_CHECK_EQUAL(bunch_train.get_bunches().at(2), bunches.at(2));
}

BOOST_FIXTURE_TEST_CASE(get_spacings1, Bunches_fixture)
{
    std::vector<double > separations;
    const double bunch_separation1 = 1.7;
    separations.push_back(bunch_separation1);
    const double bunch_separation2 = 3.4;
    separations.push_back(bunch_separation2);
    Bunch_train bunch_train(bunches, separations);

    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(0), bunch_separation1);
    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(1), bunch_separation2);
}

BOOST_FIXTURE_TEST_CASE(get_spacings2, Bunches_fixture)
{
    const double bunch_separation = 1.7;
    Bunch_train bunch_train(bunches, bunch_separation);

    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(0), bunch_separation);
    BOOST_CHECK_EQUAL(bunch_train.get_spacings().at(1), bunch_separation);
}
