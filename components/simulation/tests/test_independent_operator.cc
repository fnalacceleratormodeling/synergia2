#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/operator.h"
#include "components/simulation/lattice_simulator.h"
#include "bunch_fixture.h"
#include "lattice_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;
const int map_order = 2;

BOOST_FIXTURE_TEST_CASE(construct, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Independent_operator independent_operator("test",
            lattice_simulator.get_operation_extractor_map_sptr());
}

BOOST_FIXTURE_TEST_CASE(append_slice, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Independent_operator independent_operator("test",
            lattice_simulator.get_operation_extractor_map_sptr());

    Lattice_element_sptr element_sptr = lattice_sptr->get_elements().front();
    double length = element_sptr->get_length();
    Lattice_element_slice_sptr first_half(new Lattice_element_slice(
            *element_sptr, 0.0, 0.5 * length));
    independent_operator.append_slice(first_half);
    Lattice_element_slice_sptr second_half(new Lattice_element_slice(
            *element_sptr, 0.5 * length, length));
    independent_operator.append_slice(second_half);
}

BOOST_FIXTURE_TEST_CASE(get_slices, Lattice_fixture)
{
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Independent_operator independent_operator("test",
            lattice_simulator.get_operation_extractor_map_sptr());

    Lattice_element_sptr element_sptr = lattice_sptr->get_elements().front();
    double length = element_sptr->get_length();
    Lattice_element_slice_sptr first_half(new Lattice_element_slice(
            *element_sptr, 0.0, 0.5 * length));
    independent_operator.append_slice(first_half);
    Lattice_element_slice_sptr second_half(new Lattice_element_slice(
            *element_sptr, 0.5 * length, length));
    independent_operator.append_slice(second_half);

    Lattice_element_slices slices(independent_operator.get_slices());
    BOOST_CHECK_EQUAL(slices.size(), 2);
}

BOOST_FIXTURE_TEST_CASE(apply, Bunch_fixture)
{
    Lattice_fixture l;
    Lattice_simulator lattice_simulator(l.lattice_sptr, map_order);
    Independent_operator independent_operator("test",
            lattice_simulator.get_operation_extractor_map_sptr());

    Lattice_element_sptr element_sptr = l.lattice_sptr->get_elements().front();
    double length = element_sptr->get_length();
    Lattice_element_slice_sptr first_half(new Lattice_element_slice(
            *element_sptr, 0.0, 0.5 * length));
    independent_operator.append_slice(first_half);
    Lattice_element_slice_sptr second_half(new Lattice_element_slice(
            *element_sptr, 0.5 * length, length));
    independent_operator.append_slice(second_half);
    lattice_simulator.construct_sliced_chef_beamline(
            independent_operator.get_slices());

    Operators operators;
    Collective_operator_sptr dummy1(new Collective_operator("dummy1"));
    operators.push_back(dummy1);
    Collective_operator_sptr dummy2(new Collective_operator("dummy2"));
    operators.push_back(dummy2);
    Collective_operator_sptr dummy3(new Collective_operator("dummy3"));
    operators.push_back(dummy3);

    independent_operator.apply(bunch, operators);
}
