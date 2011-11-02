#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/lattice_simulator.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
}

BOOST_FIXTURE_TEST_CASE(construct_sliced_chef_beamline, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*(*it), 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*(*it), 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.construct_sliced_chef_beamline(slices);
}

BOOST_FIXTURE_TEST_CASE(get_map_order, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    BOOST_CHECK_EQUAL(lattice_simulator.get_map_order(), map_order);
}

BOOST_FIXTURE_TEST_CASE(get_operation_extractor_map_sptr, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    std::list<std::string >
            names(
                    lattice_simulator.get_operation_extractor_map_sptr()->get_extractor_names());
    std::list<std::string > expected_names;
    expected_names.push_back(default_operation_extractor_name);
    expected_names.push_back(chef_mixed_operation_extractor_name);
    expected_names.push_back(chef_map_operation_extractor_name);
    expected_names.push_back(chef_propagate_operation_extractor_name);

    BOOST_CHECK_EQUAL(names.size(), expected_names.size());
    names.sort();
    expected_names.sort();
    for (std::list<std::string >::iterator it = names.begin(), expected_it =
            expected_names.begin(); it != names.end(); ++it, ++expected_it) {
        BOOST_CHECK_EQUAL((*it), (*expected_it));
    }

    BOOST_CHECK_EQUAL(lattice_simulator.get_map_order(), map_order);
}

BOOST_FIXTURE_TEST_CASE(get_lattice_sptr, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.get_lattice_sptr();
}

BOOST_FIXTURE_TEST_CASE(get_chef_lattice_sptr, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.get_chef_lattice_sptr();
}

BOOST_FIXTURE_TEST_CASE(calculate_element_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    JetParticle::createStandardEnvironments(map_order);
    lattice_simulator.calculate_element_lattice_functions();
}

BOOST_FIXTURE_TEST_CASE(calculate_slice_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    JetParticle::createStandardEnvironments(map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*(*it), 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*(*it), 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.construct_sliced_chef_beamline(slices);
    lattice_simulator.calculate_slice_lattice_functions();
}

BOOST_FIXTURE_TEST_CASE(get_element_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    JetParticle::createStandardEnvironments(map_order);
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        lattice_simulator.get_lattice_functions(*(*it));
    }
}

BOOST_FIXTURE_TEST_CASE(get_slice_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    JetParticle::createStandardEnvironments(map_order);
    Lattice_element_slices slices;
    for (Lattice_elements::const_iterator it =
            lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        double length = (*it)->get_length();
        Lattice_element_slice_sptr first_half(
                new Lattice_element_slice(*(*it), 0.0, 0.5 * length));
        Lattice_element_slice_sptr second_half(
                new Lattice_element_slice(*(*it), 0.5 * length, length));
        slices.push_back(first_half);
        slices.push_back(second_half);
    }
    lattice_simulator.construct_sliced_chef_beamline(slices);
    for (Lattice_element_slices::iterator it = slices.begin(); it
            != slices.end(); ++it) {
        lattice_simulator.get_lattice_functions(*(*it));
    }
}

BOOST_FIXTURE_TEST_CASE(get_tunes, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    JetParticle::createStandardEnvironments(map_order);
    const double tolerance = 1.0e-3;
    const double expected_horizontal_tune = 0.70859;
    const double expected_vertical_tune = 0.00865009;
    BOOST_CHECK_CLOSE(lattice_simulator.get_horizontal_tune(),
            expected_horizontal_tune, tolerance);
    BOOST_CHECK_CLOSE(lattice_simulator.get_vertical_tune(),
            expected_vertical_tune, tolerance);
}

BOOST_FIXTURE_TEST_CASE(adjust_tunes, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    JetParticle::createStandardEnvironments(map_order);

    Lattice_elements horizontal_correctors, vertical_correctors;
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "quadrupole") {
            if ((*it)->get_double_attribute("k1") > 0.0) {
                horizontal_correctors.push_back(*it);
            } else {
                vertical_correctors.push_back(*it);
            }
        }
    }
    const double new_horizontal_tune = 0.69;
    const double new_vertical_tune = 0.15;
    lattice_simulator.adjust_tunes(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors);
   //jfa: need to check new tunes
}
