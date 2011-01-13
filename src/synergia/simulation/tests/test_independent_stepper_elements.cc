#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture2)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int steps_per_element = 1;
    Independent_stepper_elements stepper(lattice_simulator, steps_per_element);
    BOOST_CHECK_EQUAL(stepper.get_steps().size(),
            lattice_sptr->get_elements().size());
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture2)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int steps_per_element = 0;
    bool caught_error = false;
    try {
        Independent_stepper_elements stepper(lattice_simulator,
                steps_per_element);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture2)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int steps_per_element = 2;
    Independent_stepper_elements stepper(lattice_simulator, steps_per_element);
    //assume Lattice_fixture2 has a single zero-length element
    BOOST_CHECK_EQUAL(stepper.get_steps().size(),
            (lattice_sptr->get_elements().size()-1)*steps_per_element + 1);
}

BOOST_FIXTURE_TEST_CASE(construct17, Lattice_fixture2)
{
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    const int steps_per_element = 17;
    Independent_stepper_elements stepper(lattice_simulator, steps_per_element);
    //assume Lattice_fixture2 has a single zero-length element
    BOOST_CHECK_EQUAL(stepper.get_steps().size(),
            (lattice_sptr->get_elements().size()-1)*steps_per_element + 1);
}
