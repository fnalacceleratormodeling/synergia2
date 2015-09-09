#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/hdf5_file.h"
#include "linear_cylindrical_charge_density.h"
#include "space_charge_bunch_fixtures.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

Distributed_rectangular_grid_sptr
get_linear_cylindrical_rho2(Space_charge_2d_open_hockney & space_charge,
        Bunch & bunch, double r0, double z_period)
{
    // This is a roundabout way to set rho. We just duplicate the
    // get_global_charge_density2 test and change the values afterward
    Rectangular_grid_sptr local_rho = space_charge.get_local_charge_density(
            bunch); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2(*local_rho, bunch.get_comm_sptr()); // [C/m^3]
    std::vector<int > doubled_shape(rho2->get_domain().get_grid_shape());
    for (int i = 0; i < doubled_shape[0]; ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            rho2->get_grid_points_2dc()[i][j] = 0.0;
        }
    }
    for (int k = 0; k < doubled_shape[2]; ++k) {
        rho2->get_grid_points_1d()[k] = 0.0;
    }
    std::vector<int > nondoubled_shape(
            space_charge.get_domain().get_grid_shape());
    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<double > h(rho2->get_domain().get_cell_size());
    double factor = Q / bunch.get_total_num() / (h[0] * h[1] * h[2]);
    for (int i = 0; i < doubled_shape[0]; ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            for (int k = 0; k < doubled_shape[2]; ++k) {
                double x, y, z;
                local_rho->get_domain().get_cell_coordinates(i, j, k, x,
                        y, z);
                double r = std::sqrt(x * x + y * y);
                double val = linear_cylindrical_charge_density(lambda, r, r0);
                rho2->get_grid_points_2dc()[i][j] += val;
                rho2->get_grid_points_1d()[k] += val / factor;
            }
        }
    }
    return rho2;
}

BOOST_FIXTURE_TEST_CASE(get_local_force2_exact_rho, Cylindrical_bunch_fixture)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    Space_charge_2d_open_hockney space_charge(comm_sptr, grid_shape, true,
            false, z_period, true);
    Distributed_rectangular_grid_sptr rho2(get_linear_cylindrical_rho2(
            space_charge, bunch, r0, z_period));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    double q = bunch.get_particle_charge() * pconstants::e;
    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain().get_grid_shape());
    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    for (int component = 0; component < 2; ++component) {
        #if USE_DEBUG
        std::string filename;
        if (component == 0) {
            filename = "fx.h5";
        } else if (component == 1) {
            filename = "fy.h5";
        }
        Hdf5_file f(filename);
        f.write(local_force2->get_grid_points(), "fn");
        f.write(local_force2->get_normalization(), "fnnorm");
        #endif
        for (int i = std::max(nondoubled_shape[0] / 2,
                local_force2->get_lower()); i < std::min(
                        3 * nondoubled_shape[0] / 2, local_force2->get_upper());
                ++i) {
            for (int j = nondoubled_shape[1] / 2; j < 3 * nondoubled_shape[1]
                    / 2; ++j) {
//                for (int k = 0; k < doubled_shape[2]; ++k) {
//                    int i = nondoubled_shape[0];
//                    int j = nondoubled_shape[1];
                    int k = nondoubled_shape[2] / 2;
                    double x, y, z;
                    local_force2->get_domain().get_cell_coordinates(i,
                            j, k, x, y, z);
                    double r = std::sqrt(x * x + y * y);
                    double var = 0.0;
                    double force2_exact_ijk, force2_calc_ijk = 0.0;
                    if (component == 0) {
                        var = x;
                        force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].real()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    } else if (component == 1) {
                        var = y;
                        force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].imag()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    }
                    force2_exact_ijk
                            = linear_cylindrical_electric_force_component(q,
                                    lambda, r, r0, var) ;
//                    const double tiny = 1.0e-8;
                    double fractional_error;
                    //if (std::abs(force2_exact_ijk) < tiny) {
                    //    fractional_error = force2_calc_ijk - force2_exact_ijk;
                    //} else {
                        fractional_error = (force2_calc_ijk - force2_exact_ijk)
                                / force2_exact_ijk;
                    //}
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                    #if PRINT_FORCE
                    std::cout << x << "  " << y << "  " << z << "  "
                            << force2_exact_ijk << "  "
                            << force2_calc_ijk << "  "
                            << fractional_error << std::endl;
                    #endif
                    #if PRINT_ELECTRIC_FIELD
                    std::cout << x << "  " << y << "  " << z << "  "
                            << force2_exact_ijk / q << "  "
                            << force2_calc_ijk / q << "  "
                            << fractional_error << std::endl;
                    #endif
//                }
            }
        }
        std::cout << "get_local_force2_exact_rho" << std::endl;
        std::cout << "max_fractional_error = " << max_fractional_error
                << std::endl;
        std::cout << "min_fractional_error = " << min_fractional_error
                << std::endl;
        // on the development machine, I get
        //        max_fractional_error = 0.0183388
        //        min_fractional_error = -0.0705193
        //        max_fractional_error = 0.0183388
        //        min_fractional_error = -0.0705193
        const double field_tolerance[] = { 0.1, 0.1 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}

#if USE_POPULATE_LINEAR_CYLINDER
BOOST_FIXTURE_TEST_CASE(get_local_force2_particles, Cylindrical_bunch_fixture)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    double stdqp = 1.0e-10; // not important...
    populate_linear_cylinder(distribution, bunch, r0, z_period, stdqp, stdqp,
            stdqp);
    Space_charge_2d_open_hockney space_charge(comm_sptr, grid_shape, true,
            false, z_period, true);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch));
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho, comm_sptr));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain().get_grid_shape());
    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    for (int component = 0; component < 2; ++component) {
        #if USE_DEBUG
        std::string filename;
        if (component == 0) {
            filename = "fx.h5";
        } else if (component == 1) {
            filename = "fy.h5";
        }
        Hdf5_file f(filename);
        f.write(local_force2->get_grid_points(), "fn");
        f.write(local_force2->get_normalization(), "fnnorm");
        #endif
        for (int i = std::max(nondoubled_shape[0] / 2,
                local_force2->get_lower()); i < std::min(
                        3 * nondoubled_shape[0] / 2, local_force2->get_upper());
                ++i) {
            for (int j = nondoubled_shape[1] / 2; j < 3 * nondoubled_shape[1]
                    / 2; ++j) {
//                for (int k = 0; k < doubled_shape[2]; ++k) {
//                    int i = nondoubled_shape[0];
//                    int j = nondoubled_shape[1];
                    int k = nondoubled_shape[2] / 2;
                    double x, y, z;
                    local_force2->get_domain().get_cell_coordinates(i,
                            j, k, x, y, z);
                    double r = std::sqrt(x * x + y * y);
                    double var;
                    double force2_exact_ijk, force2_calc_ijk;
                    if (component == 0) {
                        var = x;
                        force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].real()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    } else if (component == 1) {
                        var = y;
                        force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].imag()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    }
                    force2_exact_ijk
                            = linear_cylindrical_electric_force_component(q,
                                    lambda, r, r0, var);
                    const double tiny = 1.0e-8;
                    double fractional_error;
                    if (std::abs(force2_exact_ijk) < tiny) {
                        fractional_error = force2_calc_ijk - force2_exact_ijk;
                    } else {
                        fractional_error = (force2_calc_ijk - force2_exact_ijk)
                                / force2_exact_ijk;
                    }
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                    #if PRINT_FORCE
                    std::cout << x << "  " << y << "  " << z << "  "
                            << force2_exact_ijk << "  "
                            << force2_calc_ijk << "  "
                            << fractional_error << std::endl;
                    #endif
                    #if PRINT_ELECTRIC_FIELD
                    std::cout << x << "  " << y << "  " << z << "  "
                            << force2_exact_ijk / q << "  "
                            << force2_calc_ijk / q << "  "
                            << fractional_error << std::endl;
                    #endif
//                }
            }
        }
        std::cout << "get_local_force2_particles" << std::endl;
        std::cout << "max_fractional_error = " << max_fractional_error
                << std::endl;
        std::cout << "min_fractional_error = " << min_fractional_error
                << std::endl;
        // on the development machine, I get
        //    max_fractional_error = 0.0779591
        //    min_fractional_error = -0.0980308
        //    max_fractional_error = 0.0779591
        //    min_fractional_error = -0.0980308
        //    max_fractional_error = 25655.1
        //    min_fractional_error = -25655.1
        const double field_tolerance[] = { 1.0e-1, 1.0e-1 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}
#endif

BOOST_FIXTURE_TEST_CASE(get_global_force2_exact_rho, Cylindrical_bunch_fixture)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    Space_charge_2d_open_hockney space_charge(comm_sptr, grid_shape, true, 
            false, z_period, true);
    Distributed_rectangular_grid_sptr rho2(get_linear_cylindrical_rho2(
            space_charge, bunch, r0, z_period));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    for (int component = 0; component < 2; ++component) {
        Rectangular_grid_sptr force2(
                space_charge.get_global_electric_force2(*local_force2)); // [N]
        for (int i = local_force2->get_lower(); i < local_force2->get_upper();
                ++i) {
            for (int j = 0; j
                    < local_force2->get_domain().get_grid_shape()[1];
                    ++j) {
                // Fx
                BOOST_CHECK_CLOSE(
                        local_force2->get_grid_points_2dc()[i][j].real(),
                        force2->get_grid_points_2dc()[i][j].real(),
                        tolerance);
                // Fy
                BOOST_CHECK_CLOSE(
                        local_force2->get_grid_points_2dc()[i][j].imag(),
                        force2->get_grid_points_2dc()[i][j].imag(),
                        tolerance);
            }
        }
    }
}
