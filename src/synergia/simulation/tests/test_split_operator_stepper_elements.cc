#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const int map_order = 1;
const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 1;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    BOOST_CHECK_EQUAL(stepper.get_steps().size(),
            lattice_sptr->get_elements().size());
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 0;
    bool caught_error = false;
    try {
        Split_operator_stepper_elements stepper(lattice_simulator,
                space_charge, steps_per_element);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 2;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    //assume Lattice_fixture2 has a single zero-length element
    BOOST_CHECK_EQUAL(stepper.get_steps().size(),
            (lattice_sptr->get_elements().size()-1)*steps_per_element + 1);
}

BOOST_FIXTURE_TEST_CASE(construct17, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 17;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    //assume Lattice_fixture2 has a single zero-length element
    BOOST_CHECK_EQUAL(stepper.get_steps().size(),
            (lattice_sptr->get_elements().size()-1)*steps_per_element + 1);
}

void
verify_steps(Split_operator_stepper_elements & stepper, int steps_per_element)
{
    // jfa: remove print command after debugging
    stepper.print();
    const double end = -1;
    double last_right = end;
    int slices_per_element = 2 * steps_per_element;
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
                Lattice_element_slices slices(boost::static_pointer_cast<
                        Independent_operator >(*oit)->get_slices());
                // Test 4: only one slice per operator
                BOOST_CHECK(slices.size() == 1);
                for (Lattice_element_slices::iterator slit = slices.begin(); slit
                        != slices.end(); ++slit) {
                    double slice_length = (*slit)->get_right()
                            - (*slit)->get_left();
                    // Test 5: slices should have the correct length
                    const double slice_length_tolerance = 1.0e-12;
                    BOOST_CHECK_CLOSE(slice_length,
                            (*slit)->get_lattice_element().get_length() / slices_per_element,
                            slice_length_tolerance);
                    // Test 6: slices should be continuous
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
            }
        }
    }
    // Test 7: we should end at the end of an element
    BOOST_CHECK(last_right == end);
}

BOOST_FIXTURE_TEST_CASE(verify_steps1, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 1;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    verify_steps(stepper, steps_per_element);
}

BOOST_FIXTURE_TEST_CASE(verify_steps2, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 2;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    verify_steps(stepper, steps_per_element);
}

BOOST_FIXTURE_TEST_CASE(verify_steps17, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 17;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    verify_steps(stepper, steps_per_element);
}

BOOST_FIXTURE_TEST_CASE(has_sliced_chef_beamline, Lattice_fixture2)
{
    Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int steps_per_element = 1;
    Split_operator_stepper_elements stepper(lattice_simulator, space_charge,
            steps_per_element);
    BOOST_CHECK(
            ! lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->empty());
}

