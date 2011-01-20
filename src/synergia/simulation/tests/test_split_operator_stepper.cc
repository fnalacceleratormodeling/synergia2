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
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper1(lattice_simulator, space_charge, 1);
    //    stepper1.print();
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper2(lattice_simulator, space_charge, 2);
    //    stepper2.print();
}

BOOST_FIXTURE_TEST_CASE(construct7, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper7(lattice_simulator, space_charge, 7);
    //    stepper7.print();
}

BOOST_FIXTURE_TEST_CASE(construct100, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper100(lattice_simulator, space_charge, 100);
    //    stepper100.print();
}

BOOST_FIXTURE_TEST_CASE(get_steps, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper100(lattice_simulator, space_charge, 100);

    BOOST_CHECK_EQUAL(stepper100.get_steps().size(), 100);
}

BOOST_FIXTURE_TEST_CASE(has_sliced_chef_beamline, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper stepper2(lattice_simulator, space_charge, 2);

    BOOST_CHECK(
            ! lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->empty());
}
