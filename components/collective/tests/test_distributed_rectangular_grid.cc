#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/collective/distributed_rectangular_grid.h"
#include "rectangular_grid_domain_fixture.h"

const double tolerance = 1.0e-12;
int grid_midpoint0 = grid_size0 / 2;

BOOST_FIXTURE_TEST_CASE(construct1, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(physical_size,
            physical_offset, grid_shape, is_periodic, 0, grid_size0);
}

BOOST_FIXTURE_TEST_CASE(construct2, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(
            rectangular_grid_domain_sptr, 0, grid_size0);
}

BOOST_FIXTURE_TEST_CASE(get_domain_sptr, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(
            rectangular_grid_domain_sptr, 0, grid_size0);
    BOOST_CHECK_EQUAL(rectangular_grid_domain_sptr,
            distributed_rectangular_grid.get_domain_sptr());
}

BOOST_FIXTURE_TEST_CASE(periodic_true, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(physical_size,
            physical_offset, grid_shape, true, 0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid.get_domain_sptr()->is_periodic(), true);
}

BOOST_FIXTURE_TEST_CASE(periodic_false, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(physical_size,
            physical_offset, grid_shape, false, 0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid.get_domain_sptr()->is_periodic(), false);
}

BOOST_FIXTURE_TEST_CASE(get_lower, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid1(physical_size,
            physical_offset, grid_shape, false, 0, grid_midpoint0);
    Distributed_rectangular_grid distributed_rectangular_grid2(physical_size,
            physical_offset, grid_shape, false, grid_midpoint0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid1.get_lower(), 0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid2.get_lower(), grid_midpoint0);
}

BOOST_FIXTURE_TEST_CASE(get_upper, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid1(physical_size,
            physical_offset, grid_shape, false, 0, grid_midpoint0);
    Distributed_rectangular_grid distributed_rectangular_grid2(physical_size,
            physical_offset, grid_shape, false, grid_midpoint0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid1.get_upper(), grid_midpoint0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid2.get_upper(), grid_size0);
}

BOOST_FIXTURE_TEST_CASE(get_lower_guard, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid1(physical_size,
            physical_offset, grid_shape, false, 0, grid_midpoint0);
    Distributed_rectangular_grid distributed_rectangular_grid2(physical_size,
            physical_offset, grid_shape, false, grid_midpoint0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid1.get_lower_guard(), 0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid2.get_lower_guard(), grid_midpoint0 - 1);
}

BOOST_FIXTURE_TEST_CASE(get_upper_guard, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid1(physical_size,
            physical_offset, grid_shape, false, 0, grid_midpoint0);
    Distributed_rectangular_grid distributed_rectangular_grid2(physical_size,
            physical_offset, grid_shape, false, grid_midpoint0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid1.get_upper_guard(), grid_midpoint0 + 1);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid2.get_upper_guard(), grid_size0);
}

BOOST_FIXTURE_TEST_CASE(get_lower_guard_periodic, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid1(physical_size,
            physical_offset, grid_shape, true, 0, grid_midpoint0);
    Distributed_rectangular_grid distributed_rectangular_grid2(physical_size,
            physical_offset, grid_shape, true, grid_midpoint0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid1.get_lower_guard(), -1);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid2.get_lower_guard(), grid_midpoint0 - 1);
}

BOOST_FIXTURE_TEST_CASE(get_upper_guard_periodic, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid1(physical_size,
            physical_offset, grid_shape, true, 0, grid_midpoint0);
    Distributed_rectangular_grid distributed_rectangular_grid2(physical_size,
            physical_offset, grid_shape, true, grid_midpoint0, grid_size0);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid1.get_upper_guard(), grid_midpoint0 + 1);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid2.get_upper_guard(), grid_size0 + 1);
}

BOOST_FIXTURE_TEST_CASE(get_grid_points, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(
            rectangular_grid_domain_sptr, 0, grid_size0);
    MArray3d_ref grid_points(distributed_rectangular_grid.get_grid_points());

    BOOST_CHECK_EQUAL(grid_points.shape()[0], grid_size0);
    BOOST_CHECK_EQUAL(grid_points.shape()[1], grid_size1);
    BOOST_CHECK_EQUAL(grid_points.shape()[2], grid_size2);
}

BOOST_FIXTURE_TEST_CASE(get_set_normalization, Rectangular_grid_domain_fixture)
{
    Distributed_rectangular_grid distributed_rectangular_grid(
            rectangular_grid_domain_sptr, 0, grid_size0);
    BOOST_CHECK_CLOSE(distributed_rectangular_grid.get_normalization(),
            1.0, tolerance);
    double new_norm = 123.456;
    distributed_rectangular_grid.set_normalization(new_norm);
    BOOST_CHECK_CLOSE(distributed_rectangular_grid.get_normalization(),
            new_norm, tolerance);
}

