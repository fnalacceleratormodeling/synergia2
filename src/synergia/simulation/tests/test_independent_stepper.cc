#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const int map_order = 1;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 1;
    Independent_stepper stepper(lattice_simulator, num_steps);
    //stepper.print();
}

BOOST_FIXTURE_TEST_CASE(construct_bad, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 0;
    bool caught_error = false;
    try {
        Independent_stepper stepper(lattice_simulator, num_steps);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error);
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);
}

BOOST_FIXTURE_TEST_CASE(construct17, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 17;
    Independent_stepper stepper(lattice_simulator, num_steps);
}

BOOST_FIXTURE_TEST_CASE(construct100, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Independent_stepper stepper100(lattice_simulator, 100);
    //    stepper100.print();
}

BOOST_FIXTURE_TEST_CASE(get_steps, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Independent_stepper stepper100(lattice_simulator, 100);

    BOOST_CHECK_EQUAL(stepper100.get_steps().size(), 100);
}

void
verify_steps(Independent_stepper & stepper, Lattice & lattice, int num_steps)
{
    const double end = -1;
    double last_right = end;
    double expected_step_length = lattice.get_length() / num_steps;
    for (Steps::iterator sit = stepper.get_steps().begin(); sit
            != stepper.get_steps().end(); ++sit) {
        Operators operators((*sit)->get_operators());
        double step_length = 0.0;
        // Test 1: only one operator per step
        BOOST_CHECK(operators.size() == 1);
        for (Operators::iterator oit = operators.begin(); oit
                != operators.end(); ++oit) {
            // Test 2: the operator should be an Independent_operator
            BOOST_CHECK((*oit)->get_type() == "independent");
            Lattice_element_slices slices(boost::static_pointer_cast<
                    Independent_operator >(*oit)->get_slices());
            for (Lattice_element_slices::iterator slit = slices.begin(); slit
                    != slices.end(); ++slit) {
                double slice_length = (*slit)->get_right()
                        - (*slit)->get_left();
                step_length += slice_length;
                // Test 3: slices should be continuous
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
        // Test 4: step should have the correct length
        const double step_length_tolerance = 1.0e-12;
        BOOST_CHECK_CLOSE(step_length, expected_step_length,
                step_length_tolerance);
    }
    // Test 5: we should end at the end of an element
    BOOST_CHECK(last_right == end);
}

BOOST_FIXTURE_TEST_CASE(verify_steps1, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 1;
    Independent_stepper stepper(lattice_simulator, num_steps);

    verify_steps(stepper, *lattice_sptr, num_steps);
}

BOOST_FIXTURE_TEST_CASE(verify_steps2, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);

    verify_steps(stepper, *lattice_sptr, num_steps);
}

BOOST_FIXTURE_TEST_CASE(verify_steps17, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 17;
    Independent_stepper stepper(lattice_simulator, num_steps);

    verify_steps(stepper, *lattice_sptr, num_steps);
}

BOOST_FIXTURE_TEST_CASE(has_sliced_chef_beamline, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Independent_stepper stepper2(lattice_simulator, 2);

    BOOST_CHECK(
            ! lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->empty());
}
