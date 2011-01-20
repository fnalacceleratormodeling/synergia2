#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int num_steps = 1;
    Independent_stepper stepper(lattice_simulator, num_steps);
    //stepper.print();
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int num_steps = 0;
    bool caught_error = false;
    try {
        Independent_stepper stepper(lattice_simulator,
                num_steps);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);
}

BOOST_FIXTURE_TEST_CASE(construct17, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int num_steps = 17;
    Independent_stepper stepper(lattice_simulator, num_steps);
}
