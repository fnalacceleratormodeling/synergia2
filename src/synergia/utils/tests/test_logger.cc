#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/logger.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct1)
{
    Logger logger(0);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Logger logger(0, "logger-log2");
}

BOOST_AUTO_TEST_CASE(construct3)
{
    Logger logger("logger-log3");
}

BOOST_AUTO_TEST_CASE(doit1)
{
    Logger logger(0);
    logger << "doit1" << std::endl;
    logger.flush();
}

BOOST_AUTO_TEST_CASE(doit2)
{
    Logger logger(0, std::string("logger-log2"));
    logger << "doit2" << std::endl;
    logger.flush();
}

BOOST_AUTO_TEST_CASE(doit3)
{
    Logger logger("logger-log3");
    logger << "doit3" << std::endl;
    logger.flush();
}

BOOST_AUTO_TEST_CASE(set_stream)
{
    Logger logger(0);
    logger.set_stream(std::cerr) << "should go to stderr" << std::endl;
}

