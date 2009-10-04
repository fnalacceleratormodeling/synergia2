#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/bunch.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-15;

const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 100;
const double real_num = 2.0e12;

BOOST_AUTO_TEST_CASE(construct)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
}
