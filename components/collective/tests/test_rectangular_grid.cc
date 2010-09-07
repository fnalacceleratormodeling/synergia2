#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/collective/rectangular_grid.h"
#include "rectangular_grid_domain_fixture.h"

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct1, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(physical_size, physical_offset,
            grid_shape, is_periodic);
}

BOOST_FIXTURE_TEST_CASE(construct2, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(rectangular_grid_domain_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_domain_sptr, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(rectangular_grid_domain_sptr);
    BOOST_CHECK_EQUAL(rectangular_grid_domain_sptr,
            rectangular_grid.get_domain_sptr());
}

BOOST_FIXTURE_TEST_CASE(periodic_true, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(physical_size, physical_offset,
            grid_shape, true);
    BOOST_CHECK_EQUAL(rectangular_grid.get_domain_sptr()->is_periodic(), true);
}

BOOST_FIXTURE_TEST_CASE(periodic_false, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(physical_size, physical_offset,
            grid_shape, false);
    BOOST_CHECK_EQUAL(rectangular_grid.get_domain_sptr()->is_periodic(), false);
}

BOOST_FIXTURE_TEST_CASE(get_grid_points, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(rectangular_grid_domain_sptr);
    MArray3d_ref grid_points(rectangular_grid.get_grid_points());

    BOOST_CHECK_EQUAL(grid_points.shape()[0], grid_size0);
    BOOST_CHECK_EQUAL(grid_points.shape()[1], grid_size1);
    BOOST_CHECK_EQUAL(grid_points.shape()[2], grid_size2);
}

BOOST_FIXTURE_TEST_CASE(get_set_normalization, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(physical_size, physical_offset,
                grid_shape, true);
    BOOST_CHECK_CLOSE(rectangular_grid.get_normalization(), 1.0, tolerance);
    double new_norm = 123.456;
    rectangular_grid.set_normalization(new_norm);
    BOOST_CHECK_CLOSE(rectangular_grid.get_normalization(), new_norm,
            tolerance);
}

