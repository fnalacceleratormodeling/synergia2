#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/string_utils.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const int map_order = 1;
const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 1);
}

BOOST_FIXTURE_TEST_CASE(construct_deprecated, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Split_operator_stepper stepper(lattice_simulator, space_charge, 1);
}

BOOST_FIXTURE_TEST_CASE(construct2, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 2);
}

BOOST_FIXTURE_TEST_CASE(construct7, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 7);
}

BOOST_FIXTURE_TEST_CASE(construct100, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper100(lattice_sptr, map_order, space_charge, 100);
}

BOOST_FIXTURE_TEST_CASE(get_steps, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper100(lattice_sptr, map_order, space_charge, 100);

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

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 1);
    verify_steps(stepper, *lattice_sptr, 1);
}

BOOST_FIXTURE_TEST_CASE(verify_steps2, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 2);
    verify_steps(stepper, *lattice_sptr, 2);
}

BOOST_FIXTURE_TEST_CASE(verify_steps7, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 7);
    verify_steps(stepper, *lattice_sptr, 7);
}

const double forced_diagnostics_tolerance = 1.0e-12;

void
verify_forced_diagnostics(Split_operator_stepper & stepper,
        std::string const& forced_name, Lattice_sptr lattice_sptr)
{
    double total_length = 0.0;
    double total_right_edges = 0;
    for (Steps::const_iterator it = stepper.get_steps().begin();
            it != stepper.get_steps().end(); ++it) {
        double step_length = (*it)->get_length();
        if (step_length > 0.0) {
            BOOST_CHECK(step_length > Stepper::fixed_step_tolerance);
        }
        total_length += step_length;
        int forced_element_end_count = 0;
        int first_half_count = 0;
        int second_half_count = 0;
        int zero_length_count = 0;
        bool last_right_edge = false;
        std::string last_name;
        for (Operators::iterator oit = (*it)->get_operators().begin();
                oit != (*it)->get_operators().end(); ++oit) {
            bool in_first_half = false;
            bool in_second_half = false;
            bool in_zero_length = false;
            if ((*oit)->get_name() == "first_half") {
                ++first_half_count;
                in_first_half = true;
            } else if ((*oit)->get_name() == "second_half") {
                ++second_half_count;
                in_second_half = true;
            } else if ((*oit)->get_name() == "zero_length") {
                ++zero_length_count;
                in_zero_length = true;
            }
            if (in_first_half || in_second_half || in_zero_length) {
                Lattice_element_slices slices(
                        boost::static_pointer_cast<Independent_operator >(*oit)->get_slices());
                for (Lattice_element_slices::iterator slit = slices.begin();
                        slit != slices.end(); ++slit) {
                    if ((*slit)->has_right_edge()) {
                        ++total_right_edges;
                    }
                    last_right_edge = (*slit)->has_right_edge();
                    last_name = (*slit)->get_lattice_element().get_name();
                    if ((*slit)->has_right_edge()
                            && (*slit)->get_lattice_element().has_string_attribute(
                                    Stepper::force_diagnostics_attribute)) {
                        if (!false_string(
                                (*slit)->get_lattice_element().get_string_attribute(
                                        Stepper::force_diagnostics_attribute))) {
                            forced_element_end_count += 1;
                        }
                    }
                }
            }
            if (in_first_half) {
                BOOST_CHECK(forced_element_end_count == 0);
            }
        }
        BOOST_CHECK(zero_length_count < 2);
        if (zero_length_count == 0) {
            BOOST_CHECK(first_half_count == 1);
            BOOST_CHECK(second_half_count == 1);
        }
        BOOST_CHECK(forced_element_end_count < 2);
        if (forced_element_end_count == 1) {
            BOOST_CHECK(last_right_edge);
            BOOST_CHECK(last_name == forced_name);
        }
    }
    BOOST_CHECK_CLOSE(total_length, lattice_sptr->get_length(),
            forced_diagnostics_tolerance);
    BOOST_CHECK(total_right_edges == lattice_sptr->get_elements().size());
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics0, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 2);
    verify_forced_diagnostics(stepper, "do not find me", lattice_sptr);
}


BOOST_FIXTURE_TEST_CASE(force_diagnostics1, Lattice_fixture)
{
    std::string forced_name("f");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute,
                    "true");
        }
    }
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 2);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics2, Lattice_fixture)
{
    std::string forced_name("o");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute,
                    "true");
        }
    }
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 2);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics3, Lattice_fixture)
{
    std::string forced_name("m");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute,
                    "true");
        }
    }
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 2);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics4, Lattice_fixture)
{
    std::string forced_name("f");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute,
                    "true");
        }
    }
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 7);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics5, Lattice_fixture)
{
    std::string forced_name("o");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute,
                    "true");
        }
    }
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 7);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(has_sliced_chef_beamline, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper2(lattice_sptr, map_order, space_charge, 2);

    BOOST_CHECK(
            ! stepper2.get_lattice_simulator().get_chef_lattice_sptr()->get_sliced_beamline_sptr()->empty());
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Split_operator_stepper stepper(lattice_sptr, map_order, space_charge, 1);
    xml_save(stepper, "split_operator_stepper.xml");

    Split_operator_stepper loaded;
    xml_load(loaded, "split_operator_stepper.xml");
}

BOOST_FIXTURE_TEST_CASE(serialize2_xml, Lattice_fixture)
{
    Dummy_collective_operator_sptr space_charge(
            new Dummy_collective_operator("space_charge"));

    Stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_sptr, map_order, space_charge, 1));
    xml_save(stepper_sptr, "split_operator_stepper2.xml");

    Stepper_sptr loaded;
    xml_load(loaded, "split_operator_stepper2.xml");
}
