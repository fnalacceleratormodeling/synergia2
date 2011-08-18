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
#include "gaussian_charge_density.h"
#include "space_charge_bunch_fixtures.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

Distributed_rectangular_grid_sptr
get_gaussian_rho2(Space_charge_2d_open_hockney & space_charge, Bunch & bunch,
        double sigma, double sigmaz)
{
    // This is a roundabout way to set rho. We just duplicate the
    // get_global_charge_density2 test and change the values afterward
    Rectangular_grid_sptr local_rho = space_charge.get_local_charge_density(
            bunch); // [C/m^3]

    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2(*local_rho); // [C/m^3]
    std::vector<int > doubled_shape(rho2->get_domain_sptr()->get_grid_shape());
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            rho2->get_grid_points_2dc()[i][j] = 0.0;
        }
    }
    for (int k = 0; k < doubled_shape[2]; ++k) {
        rho2->get_grid_points_1d()[k] = 0.0;
    }

    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<double > h(rho2->get_domain_sptr()->get_cell_size());
    double factor = Q / bunch.get_total_num() / (h[0] * h[1] * h[2]);
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            for (int k = 0; k < doubled_shape[2]; ++k) {
                double x, y, z;
                local_rho->get_domain_sptr()->get_cell_coordinates(i, j, k, x,
                        y, z);
                double r2 = x * x + y * y + z * z;
                //double val = gaussian_charge_density(Q, r2, sigma);
                double val = gaussian_charge_density2(Q, x, y, z, sigma, sigmaz);
                rho2->get_grid_points_2dc()[i][j] += val;
                rho2->get_grid_points_1d()[k] += val / factor;
            }
        }
    }
    return rho2;
}

#if 0
BOOST_FIXTURE_TEST_CASE(get_local_force2_exact_rho, Spherical_bunch_fixture_2d)
{
    std::vector<int > grid_shape_xyz(3);
    grid_shape_xyz[0] = grid_shape[0];
    grid_shape_xyz[1] = grid_shape[1];
    grid_shape_xyz[2] = grid_shape[2];
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma, sigmaz));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]

    // charge of the bunch
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    // charge of the single particle
    double q = bunch.get_particle_charge() * pconstants::e;
    std::vector<int > nondoubled_shape(
             space_charge.get_domain_sptr()->get_grid_shape());
    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    bool have_points = false;

    for (int component = 0; component < 2; ++component) {
        for (int i = local_force2->get_lower(); i < local_force2->get_upper();
                ++i) {
//            for (int j = 0; j < doubled_shape[1]; ++j) {
//        for (int i = nondoubled_shape[0] / 2; i < 3 * nondoubled_shape[0] / 2;
//                ++i) {
            for (int j = nondoubled_shape[1] / 2; j < 3 * nondoubled_shape[1] 
                    / 2; ++j) {
                for (int k = 0; k < doubled_shape[2]; ++k) {
                    double x, y, z;
                    space_charge.get_doubled_domain_sptr()->get_cell_coordinates(i, j, k, x, y, z);
                    double r = std::sqrt(x * x + y * y + z * z);
                    double var;
                    double local_force2_calc_ijk;
                    if (component == 0) {
                        var = x;
                        local_force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].real()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    } else if (component == 1) {
                        var = y;
                        local_force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].imag()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    }
                    double local_force2_exact_ijk 
                            = gaussian_electric_force_component(Q, q, r, sigma, 
                                    var);
                    double local_force2_exact_ijk_2
                            = gaussian_electric_force_component2(Q, q, x, y, z, 
                                    sigma, sigmaz, var);

                    double fractional_error = (local_force2_calc_ijk 
                            - local_force2_exact_ijk_2) 
                            / local_force2_exact_ijk_2;
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                    have_points = true;

                    //BOOST_CHECK_CLOSE(local_force2_calc_ijk, local_force2_exact_ijk_2, solution_tolerance);
                    //BOOST_CHECK_CLOSE(local_force2_calc_ijk, local_force2_exact_ijk, solution_tolerance);
                }
            }
        }
        if (!have_points) {
            max_fractional_error = 0.0;
            min_fractional_error = 0.0;
        }
        std::cout << "max_fractional_error = " << max_fractional_error
                << std::endl; 
        std::cout << "min_fractional_error = " << min_fractional_error
                << std::endl;
        // on the development machine, I get
        //        max_fractional_error = 0.165817
        //        min_fractional_error = -0.950593
        //        max_fractional_error = 0.165817
        //        min_fractional_error = -0.950593
        const double field_tolerance[] = { 2.0, 2.0 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}
#endif

#if 0
BOOST_FIXTURE_TEST_CASE(get_local_force2_particles, Spherical_bunch_fixture_2d)
{
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2(*local_rho); // [C/m^3]
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]

    // charge of the bunch
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    // charge of the single particle
    double q = bunch.get_particle_charge() * pconstants::e;
    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    std::vector<int > nondoubled_shape(
             space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    bool have_points = false;
    int rank = comm.get_rank();

    for (int component = 0; component < 2; ++component) {
        for (int i = local_force2->get_lower(); i < local_force2->get_upper();
                ++i) {
//            for (int j = nondoubled_shape[1] / 2; j < 3 * nondoubled_shape[1] 
//                    / 2; ++j) {
//                for (int k = 0; k < doubled_shape[2]; ++k) {
//                    int i = nondoubled_shape[0];
                    int j = nondoubled_shape[1];
                    int k = nondoubled_shape[2] / 2;
                    double x, y, z;
                    space_charge.get_doubled_domain_sptr()->get_cell_coordinates(i, j, k, x, y, z);
                    double r = std::sqrt(x * x + y * y + z * z);
                    double var;
                    double local_force2_calc_ijk;
                    if (component == 0) {
                        var = x;
                        local_force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].real()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    } else if (component == 1) {
                        var = y;
                        local_force2_calc_ijk
                                = local_force2->get_grid_points_2dc()[i][j].imag()
                                * rho2->get_grid_points_1d()[k]
                                * local_force2->get_normalization();
                    }

                    double local_force2_exact_ijk
                            = gaussian_electric_force_component(Q, q, r, sigma, 
                                    var);
                    double local_force2_exact_ijk_2
                            = gaussian_electric_force_component2(Q, q, x, y, z, 
                                    sigma, sigmaz, var);

                    double fractional_error = (local_force2_calc_ijk
                            - local_force2_exact_ijk_2)
                            / local_force2_exact_ijk_2;
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                    have_points = true;

                    #if 1 //PRINT_FORCE
                    if (rank == 0)
                    std::cout << x << "  " << y << "  " << z << "  "
                            << local_force2_exact_ijk << "  "
                            << local_force2_calc_ijk << "  "
                            << local_force2_exact_ijk_2 << "  "
                            << fractional_error << " : rank[" << rank
                            << "]" << std::endl;
                    #endif

                    //BOOST_CHECK_CLOSE(local_force2_calc_ijk, local_force2_exact_ijk_2, solution_tolerance);
                    //BOOST_CHECK_CLOSE(local_force2_calc_ijk, local_force2_exact_ijk, solution_tolerance);
//                }
//            }
        }
        if (!have_points) {
            max_fractional_error = 0.0;
            min_fractional_error = 0.0;
        }
        std::cout << "max_fractional_error = " << max_fractional_error
                << std::endl; 
        std::cout << "min_fractional_error = " << min_fractional_error
                << std::endl;
        // on the development machine, I get (on one run)
        //        max_fractional_error = 0.198425
        //        min_fractional_error = -0.947536
        //        max_fractional_error = 0.237864
        //        min_fractional_error = -0.947536
        const double field_tolerance[] = { 2.0, 2.0 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}
#endif

BOOST_FIXTURE_TEST_CASE(get_global_force2_exact_rho, Spherical_bunch_fixture_2d)
{
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma, sigmaz));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    Rectangular_grid_sptr global_force2(
            space_charge.get_global_electric_force2(*local_force2));

    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    for (int i = local_force2->get_lower(); i < local_force2->get_upper(); 
            ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            // Fx
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].real(), 
                    global_force2->get_grid_points_2dc()[i][j].real(), 
                    tolerance);
            // Fy
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].imag(),
                    global_force2->get_grid_points_2dc()[i][j].imag(),
                    tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_force2_gatherv_bcast_exact_rho, 
        Spherical_bunch_fixture_2d)
{
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma, sigmaz));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    Rectangular_grid_sptr global_force2(
            space_charge.get_global_electric_force2_gatherv_bcast(
                    *local_force2));

    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    for (int i = local_force2->get_lower(); i < local_force2->get_upper();
            ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            // Fx
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].real(),    
                    global_force2->get_grid_points_2dc()[i][j].real(),    
                    tolerance);
            // Fy
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].imag(),
                    global_force2->get_grid_points_2dc()[i][j].imag(),
                    tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_force2_allgatherv_exact_rho, 
        Spherical_bunch_fixture_2d)
{
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma, sigmaz));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    Rectangular_grid_sptr global_force2(
            space_charge.get_global_electric_force2_allgatherv(
                    *local_force2));

    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());

    for (int i = local_force2->get_lower(); i < local_force2->get_upper();
            ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            // Fx
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].real(),    
                    global_force2->get_grid_points_2dc()[i][j].real(),    
                    tolerance);
            // Fy
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].imag(),
                    global_force2->get_grid_points_2dc()[i][j].imag(),
                    tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_force2_allreduce_exact_rho,
        Spherical_bunch_fixture_2d)
{
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma, sigmaz));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    Rectangular_grid_sptr global_force2(
            space_charge.get_global_electric_force2_allreduce(*local_force2));

    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    for (int i = local_force2->get_lower(); i < local_force2->get_upper();
            ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            // Fx
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].real(),
                    global_force2->get_grid_points_2dc()[i][j].real(),
                    tolerance);
            // Fy
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].imag(),
                    global_force2->get_grid_points_2dc()[i][j].imag(),
                    tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_force2, Spherical_bunch_fixture_2d)
{
    Space_charge_2d_open_hockney space_charge(comm, grid_shape);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2(*local_rho); // [C/m^3]
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr local_force2(
            space_charge.get_local_force2(*rho2, *G2)); // [N]
    Rectangular_grid_sptr global_force2(
            space_charge.get_global_electric_force2(*local_force2));

    std::vector<int > doubled_shape(
             space_charge.get_doubled_domain_sptr()->get_grid_shape());
    for (int i = local_force2->get_lower(); i < local_force2->get_upper();
            ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            // Fx
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].real(),
                    global_force2->get_grid_points_2dc()[i][j].real(),
                    tolerance);
            // Fy
            BOOST_CHECK_CLOSE(
                    local_force2->get_grid_points_2dc()[i][j].imag(),
                    global_force2->get_grid_points_2dc()[i][j].imag(),
                    tolerance);
        }
    }
}
