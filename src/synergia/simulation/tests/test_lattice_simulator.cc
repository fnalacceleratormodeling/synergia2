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

BOOST_FIXTURE_TEST_CASE(set_slices, Lattice_fixture)
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
    lattice_simulator.set_slices(slices);
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

BOOST_AUTO_TEST_CASE(update)
{
    const double quad_length = 0.2;
    const double quad_strength = 3.2;
    const double drift_length = 3.0;
    const double bend_length = 4.0;

    Lattice_element f("quadrupole", "f");
    f.set_double_attribute("l", quad_length);
    f.set_double_attribute("k1", quad_strength);
    Lattice_element o("drift", "o");
    o.set_double_attribute("l", drift_length);
    Lattice_element d("quadrupole", "d");
    d.set_double_attribute("l", quad_length);
    d.set_double_attribute("k1", quad_strength);

    Lattice_sptr lattice_sptr(new Lattice(name));
    lattice_sptr->append(f);
    lattice_sptr->append(o);
    lattice_sptr->append(d);
    lattice_sptr->append(o);

    const int charge = pconstants::proton_charge;
    const double mass = pconstants::mp;
    const double total_energy = 125.0;
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(charge, four_momentum);
    lattice_sptr->set_reference_particle(reference_particle);

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
    lattice_simulator.set_slices(slices);

    double orig_quad_strength;
    for (beamline::deep_iterator
            it =
                    lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_begin(); it
            != lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_end(); ++it) {
        if (std::string((*it)->Type()) == "quadrupole") {
            orig_quad_strength = (*it)->Strength();
        }
    }

    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "quadrupole") {
            (*it)->set_double_attribute("k1", 2 * quad_strength);
        }
    }

    lattice_simulator.update();

    double new_quad_strength;
    for (beamline::deep_iterator
            it =
                    lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_begin(); it
            != lattice_simulator.get_chef_lattice_sptr()->get_sliced_beamline_sptr()->deep_end(); ++it) {
        if (std::string((*it)->Type()) == "quadrupole") {
            new_quad_strength = (*it)->Strength();
        }
    }
    BOOST_CHECK_CLOSE(new_quad_strength, 2*orig_quad_strength, tolerance);
}

BOOST_FIXTURE_TEST_CASE(calculate_element_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    lattice_simulator.calculate_element_lattice_functions();
}

BOOST_FIXTURE_TEST_CASE(calculate_slice_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
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
    lattice_simulator.set_slices(slices);
    lattice_simulator.calculate_slice_lattice_functions();
}

BOOST_FIXTURE_TEST_CASE(get_element_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin(); it
            != lattice_sptr->get_elements().end(); ++it) {
        lattice_simulator.get_lattice_functions(*(*it));
    }
}

BOOST_FIXTURE_TEST_CASE(get_slice_lattice_functions, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
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
    lattice_simulator.set_slices(slices);
    for (Lattice_element_slices::iterator it = slices.begin(); it
            != slices.end(); ++it) {
        lattice_simulator.get_lattice_functions(*(*it));
    }
}

BOOST_FIXTURE_TEST_CASE(get_tunes, Fobodobo_sbend_fixture)
{
    const int map_order = 1;
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
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
    const double tolerance = 1.0e-6;
    lattice_simulator.adjust_tunes(new_horizontal_tune, new_vertical_tune,
            horizontal_correctors, vertical_correctors, tolerance);
    BOOST_CHECK(std::abs(lattice_simulator.get_horizontal_tune() -
                    new_horizontal_tune) < tolerance);
    BOOST_CHECK(std::abs(lattice_simulator.get_vertical_tune() -
                    new_vertical_tune) < tolerance);
}

