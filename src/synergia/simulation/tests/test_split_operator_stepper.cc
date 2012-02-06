#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const int map_order = 1;
const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 1);
    //    stepper1.print();
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 2);
    //    stepper2.print();
}

BOOST_FIXTURE_TEST_CASE(construct7, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 7);
    //    stepper7.print();
}

BOOST_FIXTURE_TEST_CASE(construct100, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper100(lattice_simulator, space_charge, 100);
    //    stepper100.print();
}

BOOST_FIXTURE_TEST_CASE(get_steps, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper100(lattice_simulator, space_charge, 100);

    BOOST_CHECK_EQUAL(stepper100.get_steps().size(), 100);
}

void
verify_steps(Split_operator_stepper & stepper, Lattice & lattice, int num_steps)
{
    const double end = -1;
    double last_right = end;
    double expected_substep_length = lattice.get_length() / (2 * num_steps);
    for (Steps::iterator sit = stepper.get_steps().begin(); sit
            != stepper.get_steps().end(); ++sit) {
        Operators operators((*sit)->get_operators());
        // Test 1: three operators per step
        BOOST_CHECK(operators.size() == 3);
        int operator_index = 0;
        for (Operators::iterator oit = operators.begin(); oit
                != operators.end(); ++oit) {
            ++operator_index;
            // Test 2: operator should be of the right type
            if ((operator_index == 1) || (operator_index == 3)) {
                BOOST_CHECK((*oit)->get_type() == "independent");
            } else {
                BOOST_CHECK((*oit)->get_type() == "collective");
            }
            // Test 3: operators should have the correct names
            if (operator_index == 1) {
                BOOST_CHECK ((*oit)->get_name() == "first_half");
            } else if (operator_index == 2) {
                BOOST_CHECK ((*oit)->get_name() == "space_charge");
            } else if (operator_index == 3) {
                BOOST_CHECK ((*oit)->get_name() == "second_half");
            }
            if ((*oit)->get_type() == "independent") {
                double substep_length = 0.0;
                Lattice_element_slices
                        slices(
                                boost::static_pointer_cast<Independent_operator >(
                                        *oit)->get_slices());
                for (Lattice_element_slices::iterator slit = slices.begin(); slit
                        != slices.end(); ++slit) {
                    double slice_length = (*slit)->get_right()
                            - (*slit)->get_left();
                    substep_length += slice_length;
                    // Test 4: slices should be continuous
                    double slice_left = (*slit)->get_left();
                    if (last_right == end) {
                        BOOST_CHECK(slice_left == 0.0);
                    } else {
                        const double edge_tolerance = 1.0e-12;
                        BOOST_CHECK_CLOSE(last_right, slice_left, edge_tolerance);
                    }
                    if ((*slit)->has_right_edge()) {
                        last_right = end;
                    } else {
                        last_right = (*slit)->get_right();
                    }
                }
                // Test 5: substep should have the correct length
                const double substep_length_tolerance = 1.0e-12;
                BOOST_CHECK_CLOSE(substep_length, expected_substep_length,
                        substep_length_tolerance);
            }
        }
    }// Test 6: we should end at the end of an element
    BOOST_CHECK(last_right == end);
}

BOOST_FIXTURE_TEST_CASE(verify_steps1, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 1);
    verify_steps(stepper, *lattice_sptr, 1);
}

BOOST_FIXTURE_TEST_CASE(verify_steps2, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 2);
    verify_steps(stepper, *lattice_sptr, 2);
}

BOOST_FIXTURE_TEST_CASE(verify_steps7, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 7);
    verify_steps(stepper, *lattice_sptr, 7);
}

BOOST_FIXTURE_TEST_CASE(has_sliced_chef_beamline, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper2(lattice_simulator, space_charge, 2);

    BOOST_CHECK(
            ! lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->empty());
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 1);
    xml_save(stepper, "split_operator_stepper.xml");

    Split_operator_stepper loaded;
    xml_load(loaded, "split_operator_stepper.xml");
}

BOOST_FIXTURE_TEST_CASE(serialize2_xml, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge, 1));
    xml_save(stepper_sptr, "split_operator_stepper2.xml");

    Stepper_sptr loaded;
    xml_load(loaded, "split_operator_stepper2.xml");
}
