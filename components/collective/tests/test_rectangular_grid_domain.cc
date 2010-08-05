#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/collective/rectangular_grid_domain.h"
#include "rectangular_grid_domain_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
//BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Rectangular_grid_domain_fixture)
{
}

BOOST_FIXTURE_TEST_CASE(get_physical_size, Rectangular_grid_domain_fixture)
{
    std::vector<double > physical_size =
            rectangular_grid_domain_sptr->get_physical_size();
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(physical_size.at(i), domain_max - domain_min, tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_physical_offset, Rectangular_grid_domain_fixture)
{
    std::vector<double > physical_offset =
            rectangular_grid_domain_sptr->get_physical_offset();
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(physical_offset.at(i), domain_offset, tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_grid_shape, Rectangular_grid_domain_fixture)
{
    std::vector<int > grid_shape =
            rectangular_grid_domain_sptr->get_grid_shape();
    BOOST_CHECK_EQUAL(grid_shape.at(0), grid_size0);
    BOOST_CHECK_EQUAL(grid_shape.at(1), grid_size1);
    BOOST_CHECK_EQUAL(grid_shape.at(2), grid_size2);

}

BOOST_FIXTURE_TEST_CASE(is_periodic_, Rectangular_grid_domain_fixture)
{
    BOOST_CHECK_EQUAL(rectangular_grid_domain_sptr->is_periodic(), is_periodic);
}

BOOST_FIXTURE_TEST_CASE(get_leftmost_indices_offsets, Rectangular_grid_domain_fixture)
{
    std::vector<int > small_grid_shape(3);
    small_grid_shape[0] = 2;
    small_grid_shape[1] = 2;
    small_grid_shape[2] = 2;
    Rectangular_grid_domain rectangular_grid_domain(physical_size,
            physical_offset, small_grid_shape, is_periodic);

    int ix, iy, iz;
    double offx, offy, offz;
    double testx = domain_offset - 0.1, testy = domain_offset + 0.2, testz =
            domain_offset + 0.9;
    rectangular_grid_domain.get_leftmost_indices_offsets(testx, testy, testz,
            ix, iy, iz, offx, offy, offz);
    BOOST_CHECK_EQUAL(ix,0);
    BOOST_CHECK_CLOSE(offx, 0.4, tolerance);
    BOOST_CHECK_EQUAL(iy,0);
    BOOST_CHECK_CLOSE(offy, 0.7, tolerance);
    BOOST_CHECK_EQUAL(iz,1);
    BOOST_CHECK_CLOSE(offz, 0.4, tolerance);
}
