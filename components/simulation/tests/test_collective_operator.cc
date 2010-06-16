#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/operator.h"
#include "bunch_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Collective_operator collective_operator("test");
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Collective_operator collective_operator("test");

    Operators operators;
    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    operators.push_back(dummy1);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    operators.push_back(dummy2);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    operators.push_back(dummy3);

    collective_operator.apply(bunch, operators);
}
