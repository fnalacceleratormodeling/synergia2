#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/rectangular_grid_domain.h"
#include "rectangular_grid_domain_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include <cmath>
//BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Rectangular_grid_domain_fixture)
{
}

BOOST_AUTO_TEST_CASE(construct2)
{
    const double domain_min = -1.0;
    const double domain_max = 1.0;
    const double domain_offset = 5.0;
    const int grid_size0 = 4;
    const int grid_size1 = 5;
    const int grid_size2 = 3;

    std::vector<double > physical_size(3), physical_offset(3);
    std::vector<int > grid_shape(3);
    for (int i = 0; i < 3; ++i) {
        physical_offset[i] = domain_offset;
        physical_size[i] = domain_max - domain_min;
    }
    grid_shape[0] = grid_size0;
    grid_shape[1] = grid_size1;
    grid_shape[2] = grid_size2;
    Rectangular_grid_domain domain(physical_size, physical_offset, grid_shape);
    BOOST_CHECK(!domain.is_periodic());
}

BOOST_AUTO_TEST_CASE(construct3)
{
    const double domain_min = -1.0;
    const double domain_max = 1.0;
    const double domain_offset = 5.0;
    const int grid_size0 = 4;
    const int grid_size1 = 5;
    const int grid_size2 = 3;

    std::vector<double > physical_size(3), physical_offset(3);
    std::vector<int > grid_shape(3);
    for (int i = 0; i < 3; ++i) {
        physical_offset[i] = domain_offset;
        physical_size[i] = domain_max - domain_min;
    }
    grid_shape[0] = grid_size0;
    grid_shape[1] = grid_size1;
    grid_shape[2] = grid_size2;
    bool is_periodic(true);
    Rectangular_grid_domain domain(physical_size, grid_shape, is_periodic);
    BOOST_CHECK_EQUAL(is_periodic, domain.is_periodic());
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
   // double testx = domain_offset - 0.1, testy = domain_offset + 0.2, testz =
   //         domain_offset + 0.9;
            
    double testx = domain_offset - 0.6, testy = domain_offset + 0.2, testz =
            domain_offset + 0.9;        
    rectangular_grid_domain.get_leftmost_indices_offsets(testx, testy, testz,
            ix, iy, iz, offx, offy, offz);
  //  BOOST_CHECK_EQUAL(ix,0);
 //   BOOST_CHECK_CLOSE(offx, 0.4, tolerance);
    BOOST_CHECK_EQUAL(ix,-1);
    BOOST_CHECK_CLOSE(offx, 0.9, tolerance);
    BOOST_CHECK_EQUAL(iy,0);
    BOOST_CHECK_CLOSE(offy, 0.7, tolerance);
    BOOST_CHECK_EQUAL(iz,1);
    BOOST_CHECK_CLOSE(offz, 0.4, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_cell_coordinates, Rectangular_grid_domain_fixture)
{
    std::vector<int > small_grid_shape(3);
    small_grid_shape[0] = 2;
    small_grid_shape[1] = 3;
    small_grid_shape[2] = 4;
    physical_offset[0] = 0;
    physical_offset[1] = 0;
    physical_offset[2] = 0;
    physical_size[0] = 1.0;
    physical_size[1] = 3.4;
    physical_size[2] = 5.6;
    Rectangular_grid_domain rectangular_grid_domain(physical_size,
            physical_offset, small_grid_shape, is_periodic);

    int ix, iy, iz;
    double x, y, z;

    ix = 0;
    iy = 0;
    iz = 0;
    rectangular_grid_domain.get_cell_coordinates(ix, iy, iz, x, y, z);
    BOOST_CHECK_CLOSE(x,
            -(physical_size[0]/(1.0*small_grid_shape[0]))*0.5, tolerance);
    BOOST_CHECK_CLOSE(y,
            -(physical_size[1]/(1.0*small_grid_shape[1])), tolerance);
    BOOST_CHECK_CLOSE(z,
            -(physical_size[2]/(1.0*small_grid_shape[2]))*1.5, tolerance);

    ix = 1;
    iy = 1;
    iz = 1;
    rectangular_grid_domain.get_cell_coordinates(ix, iy, iz, x, y, z);
    BOOST_CHECK_CLOSE(x,
            (physical_size[0]/(1.0*small_grid_shape[0]))*0.5, tolerance);
    BOOST_CHECK(std::abs(y) < tolerance);
    BOOST_CHECK_CLOSE(z,
            -(physical_size[2]/(1.0*small_grid_shape[2]))*0.5, tolerance);
}
