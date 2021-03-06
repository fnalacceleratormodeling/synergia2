#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/operator.h"
#include "synergia/utils/serialization.h"
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
const double step_length = 1.23;

BOOST_AUTO_TEST_CASE(construct)
{
    Step step(step_length);
}

BOOST_AUTO_TEST_CASE(append)
{
    Step step(step_length);
    Dummy_collective_operator_sptr Dummy_collective_operator_sptr(
            new Dummy_collective_operator("test"));

    double time_fraction = 1.0;
    step.append(Dummy_collective_operator_sptr, time_fraction);
}

BOOST_AUTO_TEST_CASE(append2)
{
    Step step(step_length);
    Dummy_collective_operator dummy_collective_operator("test");

    Operators operators;
    Dummy_collective_operator_sptr dummy1(
            new Dummy_collective_operator("dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(
            new Dummy_collective_operator("dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(
            new Dummy_collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.append(operators, 1.0);
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Step step(step_length);

    Operators operators;
    Dummy_collective_operator_sptr dummy1(
            new Dummy_collective_operator("dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(
            new Dummy_collective_operator("dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(
            new Dummy_collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.append(operators, 1.0);

    const int verbosity = 3;
    Logger logger(0);
    Lattice_element o("drift", "o");
    Lattice_sptr lattice_sptr(new Lattice("test_lattice"));   
    lattice_sptr->append(o);
    Reference_particle reference_particle(1, pconstants::mp, 2.);
    lattice_sptr->set_reference_particle(reference_particle);
    Lattice_simulator lattice_simulator(lattice_sptr,1);
    Independent_stepper stepper(lattice_simulator,1);
    Diagnosticss per_operator_diagnosticss, per_operation_diagnosticss;
    step.apply(bunch, verbosity, per_operator_diagnosticss,
            per_operation_diagnosticss, stepper, logger);
}

BOOST_AUTO_TEST_CASE(get_operators)
{
    Step step(step_length);

    Operators operators;
    Dummy_collective_operator_sptr dummy1(
            new Dummy_collective_operator("dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(
            new Dummy_collective_operator("dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(
            new Dummy_collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.append(operators, 1.0);

    Operators retrieved_operators(step.get_operators());
    BOOST_CHECK_EQUAL(retrieved_operators.size(),3);
}

BOOST_AUTO_TEST_CASE(get_time_fractions)
{
    Step step(step_length);

    Dummy_collective_operator_sptr dummy1(
            new Dummy_collective_operator("dummy1"));
    step.append(dummy1, 0.5);
    Dummy_collective_operator_sptr dummy2(
            new Dummy_collective_operator("dummy2"));
    step.append(dummy2, 1.0);
    Dummy_collective_operator_sptr dummy3(
            new Dummy_collective_operator("dummy3"));
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

BOOST_AUTO_TEST_CASE(get_length)
{
    Step step(step_length);

    BOOST_CHECK_CLOSE(step.get_length(), step_length, tolerance);
}

BOOST_AUTO_TEST_CASE(print)
{
    Step step(step_length);

    Operators operators;
    Dummy_collective_operator_sptr dummy1(
            new Dummy_collective_operator("dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(
            new Dummy_collective_operator("dummy2"));
    operators.push_back(dummy2);
    Dummy_collective_operator_sptr dummy3(
            new Dummy_collective_operator("dummy3"));
    operators.push_back(dummy3);

    step.print(1);
}

BOOST_AUTO_TEST_CASE(serialize_xml)
{
    Step step(step_length);
    Operators operators;
    Dummy_collective_operator_sptr dummy1(
            new Dummy_collective_operator("dummy1"));
    operators.push_back(dummy1);
    Dummy_collective_operator_sptr dummy2(
            new Dummy_collective_operator("dummy2"));
    operators.push_back(dummy2);
    step.append(operators, 1.0);

    xml_save<Step > (step, "step.xml");

    Step loaded;
    xml_load<Step > (loaded, "step.xml");
}
