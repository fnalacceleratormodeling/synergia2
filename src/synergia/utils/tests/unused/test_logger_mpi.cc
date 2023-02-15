#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/logger.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

BOOST_AUTO_TEST_CASE(construct1)
{
    Logger logger(0);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Logger logger(0, "logger-log-mpi-2");
}

BOOST_AUTO_TEST_CASE(construct3)
{
    Logger logger("logger-log-mpi-3");
}

BOOST_AUTO_TEST_CASE(doit1)
{
    Logger logger(0);
    logger << "doit1" << std::endl;
    logger.flush();
}

BOOST_AUTO_TEST_CASE(doit11)
{
    Logger logger(1);
    logger << "doit11" << std::endl;
    logger.flush();
}

BOOST_AUTO_TEST_CASE(doit2)
{
    Logger logger(0, "logger-log-mpi-2");
    logger << "doit2" << std::endl;
    logger.flush();
}

BOOST_AUTO_TEST_CASE(doit3)
{
    Logger logger("logger-log-mpi-3");
    logger << "I am doing it 3 on rank " << Commxx().get_rank() << std::endl;
    logger.flush();
}
