#define BOOST_TEST_MAIN
#include <iostream>
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
    
    std::vector<double> cellS(3,0.);
    for (int k=0; k != 3; k++) cellS[k] = (domain_max - domain_min)/((double) a.shape()[k]);
    const double left = domain_offset - 0.5*(domain_max - domain_min);
    for (int ix = 0; ix != grid_size0; ix++) {
      const double x = left +  ix*cellS[0];
      for (int iy = 0; iy != grid_size1; iy++) {
        const double y = left +  iy*cellS[1];
        for (int iz = 0; iz != grid_size2; iz++) {
          const double z = left +  iz*cellS[2];
	  a[ix][iy][iz] = 5.0 + x*2. + y*y*3. + z*z*z*4.;
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
    for (int k=0; k !=3; k++) goodCoords[k] += 0.25*cellS[k];
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
    // 
    // Not so trivial test based on derivatives. 
    //
    const double ddx = 1.267*cellS[0];
    const double ddy = 1.453*cellS[1];
    const double y1 = domain_offset + 0.123*cellS[1];
    const double z1 = domain_offset + 0.245*cellS[2];
    const double xx1 = domain_offset - ddx;
    const double xx2 = domain_offset + ddx;
    const double dfx1 = (rectangular_grid.get_interpolated_coord(xx2, y1, z1) - 
                        rectangular_grid.get_interpolated_coord(xx1, y1, z1))/(2.0*ddx);
			
    const double dfx1Pred = 2.;
//    std::cerr << " cell size " << cellS[0] << " f2 " << rectangular_grid.get_interpolated_coord(xx2, y1, z1)
//                               << " f1 " << rectangular_grid.get_interpolated_coord(xx1, y1, z1) 
//			       << " dfx1 " << dfx1 << " Pred " << dfx1Pred << std::endl;
    BOOST_CHECK_CLOSE(dfx1, dfx1Pred,  0.1);
    //
    //			
    const double x1 = domain_offset + 0.1*cellS[0];
    const double yy1 = domain_offset - ddy;
    const double yy2 = domain_offset + ddy;
    const double dfy1 = (rectangular_grid.get_interpolated_coord(x1, yy2, z1) - 
                        rectangular_grid.get_interpolated_coord(x1, yy1, z1))/(2.0*ddy);
			
    const double dfy1Pred = 6.*domain_offset;
//    std::cerr << " cell size " << cellS[1] << " f2 " << rectangular_grid.get_interpolated_coord(x1, yy2, z1)
//                               << " f1 " << rectangular_grid.get_interpolated_coord(x1, yy1, z1) << std::endl;
    BOOST_CHECK_CLOSE(dfy1, dfy1Pred,  5.); // nonlinearities on discrete lattice. .. 
    //
    
    
}
