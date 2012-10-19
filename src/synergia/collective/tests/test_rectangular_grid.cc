#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/rectangular_grid.h"
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

BOOST_FIXTURE_TEST_CASE(get_const_domain_sptr, Rectangular_grid_domain_fixture)
{
    const Rectangular_grid rectangular_grid(rectangular_grid_domain_sptr);
    BOOST_CHECK_EQUAL(rectangular_grid_domain_sptr,
            rectangular_grid.get_domain_sptr());
}

BOOST_FIXTURE_TEST_CASE(periodic_true, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(physical_size, physical_offset,
            grid_shape, true);
    BOOST_CHECK_EQUAL(rectangular_grid.get_domain().is_periodic(), true);
}

BOOST_FIXTURE_TEST_CASE(periodic_false, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(physical_size, physical_offset,
            grid_shape, false);
    BOOST_CHECK_EQUAL(rectangular_grid.get_domain().is_periodic(), false);
}

BOOST_FIXTURE_TEST_CASE(get_const_grid_points, Rectangular_grid_domain_fixture)
{
    const Rectangular_grid rectangular_grid(rectangular_grid_domain_sptr);
    MArray3d_ref grid_points(rectangular_grid.get_grid_points());

    BOOST_CHECK_EQUAL(grid_points.shape()[0], grid_size0);
    BOOST_CHECK_EQUAL(grid_points.shape()[1], grid_size1);
    BOOST_CHECK_EQUAL(grid_points.shape()[2], grid_size2);
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

BOOST_FIXTURE_TEST_CASE(interpolate_coord, Rectangular_grid_domain_fixture)
{
    Rectangular_grid rectangular_grid(rectangular_grid_domain_sptr);
    
    MArray3d_ref a(rectangular_grid.get_grid_points());
    
    for (int ix = 0; ix != grid_size0; ix++) {
      for (int iy = 0; iy != grid_size1; iy++) {
        for (int iz = 0; iz != grid_size2; iz++) {
	  a[ix][iy][iz] = 5.0 + ix*2. + iy*iy*3. + iz*iz*iz*4.;
	}
      }
    }
    // Trivial case first, out of grid in the Y direction. 
    // 
    std::vector<double> offYCoords(3, domain_offset + 0.5*(domain_max - domain_min));
    offYCoords[1] = domain_offset + 5.0*(domain_max - domain_min);
    BOOST_CHECK_CLOSE(rectangular_grid.get_interpolated_coord(offYCoords[0], offYCoords[1], offYCoords[2]),
                                         0.0, tolerance);
    //
    // not so trivial 
    //
    std::vector<double> goodCoords(3, domain_offset);
    std::vector<double> cellS(3,0.);
    MArray3d_ref grid_points(rectangular_grid.get_grid_points());
    for (int k=0; k != 3; k++) cellS[k] = (domain_max - domain_min)/((double) grid_points.shape()[k]);
    for (int k=0; k !=3; k++) goodCoords[k] += 0.25*cellS[k];
    const double left = domain_offset - 0.5*(domain_max - domain_min);
    const double sclX = (goodCoords[0] - left) / cellS[0] - 0.5;
    const double sclY = (goodCoords[1] - left) / cellS[1] - 0.5;
    const double sclZ = (goodCoords[2] - left) / cellS[2] - 0.5;
    const int ix =  fast_int_floor(sclX);		 
    const int iy =  fast_int_floor(sclY);		 
    const int iz =  fast_int_floor(sclZ);		 
    const double offx = sclX - ix;
    const double offy = sclY - iy;
    const double offz = sclZ - iz;
    const double aVal = 
                ( (1.0 - offx) * (1.0 - offy) * (1.0 - offz) * a[ix][iy][iz]
                + (1.0 - offx) * (1.0 - offy) * offz * a[ix][iy][iz + 1] 
		+ (1.0 - offx) * offy * (1.0 - offz) * a[ix][iy + 1][iz]
		+ (1.0 - offx) * offy * offz * a[ix][iy + 1][iz + 1] 
                + offx * (1.0 - offy) * (1.0 - offz) * a[ix + 1][iy][iz] 
		+ offx * (1.0 - offy) * offz * a[ix + 1][iy][iz + 1] 
		+ offx * offy * (1.0 - offz) * a[ix + 1][iy + 1][iz] 
		+ offx * offy * offz * a[ix + 1][iy + 1][iz + 1]);
       
    BOOST_CHECK_CLOSE(rectangular_grid.get_interpolated_coord(goodCoords[0], goodCoords[1], goodCoords[2]), 
                                         aVal, tolerance);
    
}
