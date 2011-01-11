#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/step.h"
#include "synergia/simulation/operator.h"
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

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

    double time_fraction = 1.0;
    step.append(collective_operator_sptr, time_fraction);
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

    step.append(operators, 1.0);
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

    step.append(operators, 1.0);

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

    step.append(operators, 1.0);

    Operators retrieved_operators(step.get_operators());
    BOOST_CHECK_EQUAL(retrieved_operators.size(),3);
}

BOOST_AUTO_TEST_CASE(get_time_fractions)
{
    Step step;

    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    step.append(dummy1, 0.5);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    step.append(dummy2, 1.0);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    step.append(dummy3, 0.25);

    std::list<double > retrieved_time_fractions(step.get_time_fractions());
    BOOST_CHECK_EQUAL(retrieved_time_fractions.size(),3);
    std::list<double >::const_iterator it = retrieved_time_fractions.begin();
    BOOST_CHECK_CLOSE(*it, 0.5, tolerance);
    ++it;
    BOOST_CHECK_CLOSE(*it, 1.0, tolerance);
    ++it;
    BOOST_CHECK_CLOSE(*it, 0.25, tolerance);
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

