#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iostream>

#include "synergia/collective/interpolate_rectangular_zyx.h"
#include "gaussian_charge_density.h"

void
set_grid_constant(Rectangular_grid & grid, double val)
{
    for (int i = 0; i < grid.get_domain_sptr()->get_grid_shape()[0]; ++i) {
        for (int j = 0; j < grid.get_domain_sptr()->get_grid_shape()[1]; ++j) {
            for (int k = 0; k < grid.get_domain_sptr()->get_grid_shape()[2]; ++k) {
                grid.get_grid_points()[i][j][k] = val;
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(interpolate_rectangular_zyx_basic)
{
    std::vector<double > physical_size(3);
    physical_size[0] = 2.0;
    physical_size[1] = 2.0;
    physical_size[2] = 2.0;
    std::vector<double > physical_offset(3);
    physical_offset[0] = 0.0;
    physical_offset[1] = 0.0;
    physical_offset[2] = 0.0;
    std::vector<int > grid_shape(3);
    grid_shape[0] = 4;
    grid_shape[1] = 4;
    grid_shape[2] = 4;

    Rectangular_grid grid(physical_size, physical_offset, grid_shape, false);
    std::vector<double > cell_size(grid.get_domain_sptr()->get_cell_size());

    double z_left, y_left, x_left;
    grid.get_domain_sptr()->get_cell_coordinates(0, 0, 0, z_left, y_left,
            x_left);
    double z_right, y_right, x_right;
    grid.get_domain_sptr()->get_cell_coordinates(0, 0, 0, z_right, y_right,
            x_right);
    double left_val = 2.0;
    double right_val = 10.0;

    const int num_steps = 16;
    const double strict_tolerance = 1.0e-14;

    set_grid_constant(grid, 1.0);
    // outside grid
    BOOST_CHECK(interpolate_rectangular_zyx(z_left - cell_size[0],
                    y_left, x_left, grid) == 0.0);
    BOOST_CHECK(interpolate_rectangular_zyx(z_left,
                    y_left - cell_size[1], x_left, grid) == 0.0);
    BOOST_CHECK(interpolate_rectangular_zyx(z_left,
                    y_left, x_left - cell_size[2], grid) == 0.0);

    // constant-value parallel planes in third coordinate of grid
    set_grid_constant(grid, 0.0);
    grid.get_grid_points()[0][0][0] = left_val;
    grid.get_grid_points()[0][1][0] = left_val;
    grid.get_grid_points()[1][0][0] = left_val;
    grid.get_grid_points()[1][1][0] = left_val;
    grid.get_grid_points()[0][0][1] = right_val;
    grid.get_grid_points()[0][1][1] = right_val;
    grid.get_grid_points()[1][0][1] = right_val;
    grid.get_grid_points()[1][1][1] = right_val;
    double dval_ds = (right_val - left_val) / cell_size[2];
    for (double ds = 0.0; ds <= cell_size[2]; ds += cell_size[2] / num_steps) {
        double expected_val = left_val + dval_ds * ds;
        // value along edge of cube
        double val = interpolate_rectangular_zyx(x_left + ds, y_left, z_left,
                grid);
        BOOST_CHECK_CLOSE(val, expected_val, strict_tolerance);
        // value along center of cube
        val = interpolate_rectangular_zyx(x_left + ds,
                (y_left + y_right) / 2.0, (z_left + z_right) / 2.0, grid);
        BOOST_CHECK_CLOSE(val, expected_val, strict_tolerance);
    }

    // constant-value parallel planes in second coordinate of grid
    set_grid_constant(grid, 0.0);
    grid.get_grid_points()[0][0][0] = left_val;
    grid.get_grid_points()[0][0][1] = left_val;
    grid.get_grid_points()[1][0][0] = left_val;
    grid.get_grid_points()[1][0][1] = left_val;
    grid.get_grid_points()[0][1][0] = right_val;
    grid.get_grid_points()[0][1][1] = right_val;
    grid.get_grid_points()[1][1][0] = right_val;
    grid.get_grid_points()[1][1][1] = right_val;
    dval_ds = (right_val - left_val) / cell_size[1];
    for (double ds = 0.0; ds <= cell_size[1]; ds += cell_size[1] / num_steps) {
        double expected_val = left_val + dval_ds * ds;
        // value along edge of cube
        double val = interpolate_rectangular_zyx(x_left, y_left + ds, z_left,
                grid);
        BOOST_CHECK_CLOSE(val, expected_val, strict_tolerance);
        // value along center of cube
        val = interpolate_rectangular_zyx((x_left + x_right) / 2.0,
                y_left + ds, (z_left + z_right) / 2.0, grid);
        BOOST_CHECK_CLOSE(val, expected_val, strict_tolerance);
    }

    // constant-value parallel planes in first coordinate of grid
    set_grid_constant(grid, 0.0);
    grid.get_grid_points()[0][0][0] = left_val;
    grid.get_grid_points()[0][0][1] = left_val;
    grid.get_grid_points()[0][1][0] = left_val;
    grid.get_grid_points()[0][1][1] = left_val;
    grid.get_grid_points()[1][0][0] = right_val;
    grid.get_grid_points()[1][0][1] = right_val;
    grid.get_grid_points()[1][1][0] = right_val;
    grid.get_grid_points()[1][1][1] = right_val;
    dval_ds = (right_val - left_val) / cell_size[0];
    for (double ds = 0.0; ds <= cell_size[0]; ds += cell_size[0] / num_steps) {
        double expected_val = left_val + dval_ds * ds;
        // value along edge of cube
        double val = interpolate_rectangular_zyx(x_left, y_left, z_left + ds,
                grid);
        BOOST_CHECK_CLOSE(val, expected_val, strict_tolerance);
        // value along center of cube
        val = interpolate_rectangular_zyx((x_left + x_right) / 2.0, (y_left
                + y_right) / 2.0, z_left + ds, grid);
        BOOST_CHECK_CLOSE(val, expected_val, strict_tolerance);
    }

}

//BOOST_AUTO_TEST_CASE(interpolate_rectangular_zyx_gaussian)
//{
//    const double sigma = 1.7e-3;
//    const double Q = 2.3e-8;
//    std::vector<double > physical_size(3);
//    physical_size[0] = 6.0 * sigma;
//    physical_size[1] = 7.0 * sigma;
//    physical_size[2] = 8.0 * sigma;
//    std::vector<double > physical_offset(3);
//    physical_offset[0] = 0.0;
//    physical_offset[1] = 0.0;
//    physical_offset[2] = 0.0;
//    std::vector<int > grid_shape(3);
//    grid_shape[0] = 32;
//    grid_shape[1] = 41;
//    grid_shape[2] = 48;
//
//    Rectangular_grid grid(physical_size, physical_offset, grid_shape, false);
//    for (int i = 0; i < grid_shape[0]; ++i) {
//        for (int j = 0; j < grid_shape[1]; ++j) {
//            for (int k = 0; k < grid_shape[2]; ++k) {
//                double z, y, x;
//                grid.get_domain_sptr()->get_cell_coordinates(i, j, k, z, y, x);
//                double r = std::sqrt(z * z + y * y + x * x);
//                grid.get_grid_points()[i][j][k]
//                        = gaussian_electric_field_component(Q, r, sigma, x);
//            }
//        }
//    }
//    double max_fractional_error = -2.0;
//    double min_fractional_error = 2.0;
//    double limit = 2.0 * sigma;
//    double step = limit / 10.0;
////    for (double x = -limit; x < limit; x += step) {
////        for (double y = -limit; y < limit; y += step) {
//    double x = sigma;
//    double y = 0.0;
//            for (double z = -limit; z < limit; z += step) {
//                double val = interpolate_rectangular_zyx(x, y, z, grid);
//                double r = std::sqrt(z * z + y * y + x * x);
//                double eval = gaussian_electric_field_component(Q, r, sigma, z);
//                std::cout
//                    << x << ", "
//                    << y << ", "
//                    << z << ": "
//                    << val << ", "
//                    << eval << std::endl;
//                double fractional_error = (val - eval) / eval;
//                if (fractional_error > max_fractional_error) {
//                    max_fractional_error = fractional_error;
//                }
//                if (fractional_error < min_fractional_error) {
//                    min_fractional_error = fractional_error;
//                }
//
//            }
////        }
////    }
//    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
//    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;
//
//    // on the development machine, I get
//    //    const double solution_tolerance = 2.0e-2;
//    //    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
//    //    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
//
//}
