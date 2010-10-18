#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Collective_operator collective_operator("test");
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Collective_operator collective_operator("test");
    double step_length = 1.0;
    Step stub_step;

    collective_operator.apply(bunch, step_length, stub_step);
}
