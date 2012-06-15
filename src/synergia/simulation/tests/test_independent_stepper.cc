#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "lattice_fixture.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/string_utils.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const int map_order = 1;
#if 0
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
            Lattice_element_slices
                    slices(
                            boost::static_pointer_cast<Independent_operator >(
                                    *oit)->get_slices());
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
#endif

const double forced_diagnostics_tolerance = 1.0e-10;

void
verify_forced_diagnostics(Independent_stepper & stepper,
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
        int operator_count = 0;
        bool last_right_edge = false;
        std::string last_name;
        for (Operators::iterator oit = (*it)->get_operators().begin();
                oit != (*it)->get_operators().end(); ++oit) {
            ++operator_count;
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
        BOOST_CHECK(operator_count == 1);
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
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);

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
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);

    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics2, Lattice_fixture)
{
    std::string forced_name("o");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it!= lattice_sptr->get_elements().end(); ++it){
        if((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute, "true");
        }
    }
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);

    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics3, Lattice_fixture)
{
    std::string forced_name("m");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it!= lattice_sptr->get_elements().end(); ++it){
        if((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute, "true");
        }
    }
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

BOOST_FIXTURE_TEST_CASE(force_diagnostics4, Lattice_fixture)
{
    std::string forced_name("o");
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it!= lattice_sptr->get_elements().end(); ++it){
        if((*it)->get_name() == forced_name) {
            (*it)->set_string_attribute(Stepper::force_diagnostics_attribute, "true");
        }
    }
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 7;
    Independent_stepper stepper(lattice_simulator, num_steps);
    verify_forced_diagnostics(stepper, forced_name, lattice_sptr);
}

#if 0
BOOST_FIXTURE_TEST_CASE(has_sliced_chef_beamline, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    Independent_stepper stepper2(lattice_simulator, 2);

    BOOST_CHECK(
            ! lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->empty());
}

BOOST_FIXTURE_TEST_CASE(serialize_xml, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 2;
    Independent_stepper stepper(lattice_simulator, num_steps);
    xml_save(stepper, "independent_stepper.xml");

    Independent_stepper loaded;
    xml_load(loaded, "independent_stepper.xml");
}

BOOST_FIXTURE_TEST_CASE(serialize2_xml, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    const int num_steps = 2;
    Stepper_sptr stepper_sptr(
            new Independent_stepper(lattice_simulator, num_steps));
    xml_save(stepper_sptr, "independent_stepper2.xml");

    Stepper_sptr loaded;
    xml_load(loaded, "independent_stepper2.xml");
}
#endif
