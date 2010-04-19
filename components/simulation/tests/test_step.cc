#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/step.h"
#include "components/simulation/operator.h"
#include "bunch_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Step step;
}

BOOST_AUTO_TEST_CASE(append)
{
    Step step;
    Collective_operator_sptr collective_operator_sptr(new Collective_operator(
            "test"));

    step.append(collective_operator_sptr);
}

BOOST_AUTO_TEST_CASE(append2)
{
    Step step;
    Collective_operator collective_operator("test");

    Operators operators;
    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    operators.push_back(dummy1);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    operators.push_back(dummy2);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.append(operators);
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Step step;

    Operators operators;
    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    operators.push_back(dummy1);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    operators.push_back(dummy2);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.append(operators);

    step.apply(bunch);
}

BOOST_AUTO_TEST_CASE(get_operators)
{
    Step step;

    Operators operators;
    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    operators.push_back(dummy1);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    operators.push_back(dummy2);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.append(operators);

    Operators retrieved_operators(step.get_operators());
    BOOST_CHECK_EQUAL(retrieved_operators.size(),3);
}

BOOST_AUTO_TEST_CASE(print)
{
    Step step;

    Operators operators;
    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    operators.push_back(dummy1);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    operators.push_back(dummy2);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.print(1);
}

