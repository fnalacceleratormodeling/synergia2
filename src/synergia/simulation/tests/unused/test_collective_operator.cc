#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "bunch_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/serialization.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct)
{
    Dummy_collective_operator collective_operator("test");
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Dummy_collective_operator collective_operator("test");
    double step_length = 1.0;
    Step stub_step(1.0);

    const int verbosity = 4;
    Logger logger(0);
    collective_operator.apply(bunch, step_length, stub_step, verbosity, logger);
}

BOOST_AUTO_TEST_CASE(serialize_xml)
{
    Dummy_collective_operator collective_operator("test");
    xml_save<Dummy_collective_operator > (collective_operator,
            "collective_operator.xml");

    Dummy_collective_operator loaded;
    xml_load<Dummy_collective_operator > (loaded, "collective_operator.xml");

    BOOST_CHECK_EQUAL(loaded.get_name(), "test");
    BOOST_CHECK_EQUAL(loaded.get_type(), "collective");
}

BOOST_AUTO_TEST_CASE(serialize2_xml)
{
    Operator_sptr operator_sptr(new Dummy_collective_operator("test"));
    xml_save(operator_sptr, "collective_operator2.xml");

    Operator_sptr loaded;
    xml_load(loaded, "collective_operator2.xml");
}

