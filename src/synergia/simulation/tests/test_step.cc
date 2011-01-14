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
    Dummy_collective_operator_sptr Dummy_collective_operator_sptr(
            new Dummy_collective_operator("test"));

    double time_fraction = 1.0;
    step.append(Dummy_collective_operator_sptr, time_fraction);
}

BOOST_AUTO_TEST_CASE(append2)
{
    Step step;
    Dummy_collective_operator dummy_collective_operator("test");

    Operators operators;
    Dummy_collective_operator_sptr dummy1(new Dummy_collective_operator(
            "dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(new Dummy_collective_operator(
            "dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(new Dummy_collective_operator(
            "dummy3"));
    operators.push_back(dummy3);

    step.append(operators, 1.0);
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Step step;

    Operators operators;
    Dummy_collective_operator_sptr dummy1(new Dummy_collective_operator(
            "dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(new Dummy_collective_operator(
            "dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(new Dummy_collective_operator(
            "dummy3"));
    operators.push_back(dummy3);

    step.append(operators, 1.0);

    step.apply(bunch);
}

BOOST_AUTO_TEST_CASE(get_operators)
{
    Step step;

    Operators operators;
    Dummy_collective_operator_sptr dummy1(new Dummy_collective_operator(
            "dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(new Dummy_collective_operator(
            "dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(new Dummy_collective_operator(
            "dummy3"));
    operators.push_back(dummy3);

    step.append(operators, 1.0);

    Operators retrieved_operators(step.get_operators());
    BOOST_CHECK_EQUAL(retrieved_operators.size(),3);
}

BOOST_AUTO_TEST_CASE(get_time_fractions)
{
    Step step;

    Dummy_collective_operator_sptr dummy1(new Dummy_collective_operator(
            "dummy1"));
    step.append(dummy1, 0.5);
    Dummy_collective_operator_sptr dummy2(new Dummy_collective_operator(
            "dummy2"));
    step.append(dummy2, 1.0);
    Dummy_collective_operator_sptr dummy3(new Dummy_collective_operator(
            "dummy3"));
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
    Dummy_collective_operator_sptr dummy1(new Dummy_collective_operator(
            "dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(new Dummy_collective_operator(
            "dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(new Dummy_collective_operator(
            "dummy3"));
    operators.push_back(dummy3);

    step.print(1);
}

