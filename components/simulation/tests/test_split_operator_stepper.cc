#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/stepper.h"
#include "components/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Collective_operator_sptr space_charge(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    lattice_sptr->print();

    Split_operator_stepper stepper1(lattice_simulator, space_charge, 1);
    stepper1.print();

    Split_operator_stepper stepper2(lattice_simulator, space_charge, 2);
    stepper2.print();

    Split_operator_stepper stepper7(lattice_simulator, space_charge, 7);
    stepper7.print();

    Split_operator_stepper stepper10(lattice_simulator, space_charge, 10);

    stepper10.print();
}

