#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "synergia/collective/space_charge_3d_open_hockney.h"
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

// define DBGPRINT to get verbose output of test values
#define DBGPRINT 0

const double tolerance = 1.0e-12;

Distributed_rectangular_grid_sptr
get_gaussian_rho2(Space_charge_3d_open_hockney & space_charge, Bunch & bunch,
        double sigma)
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
            for (int k = 0; k < doubled_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] = 0.0;
            }
        }
    }
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    int i_min = std::max(0, rho2->get_lower());
    int i_max = std::min(nondoubled_shape[0], rho2->get_upper());
    for (int i = i_min; i < i_max; ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                local_rho->get_domain_sptr()->get_cell_coordinates(i, j, k, z,
                        y, x);
                double r2 = x * x + y * y + z * z;
                rho2->get_grid_points()[i][j][k] = gaussian_charge_density(Q,
                        r2, sigma);
            }
        }
    }
    return rho2;
}

BOOST_FIXTURE_TEST_CASE(get_scalar_field2_exact_rho, Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper(), Commxx());
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    double have_points = false;
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain_sptr()->get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y + z * z);
                double phi_exact_ijk = gaussian_electric_potential(Q, r, sigma);
                phi_exact.get_grid_points()[i][j][k] = phi_exact_ijk;
                double phi_calc_ijk = phi2->get_grid_points()[i][j][k]
                        * phi2->get_normalization();
                double fractional_error = (phi_calc_ijk - phi_exact_ijk)
                        / phi_exact_ijk;
                if (fractional_error > max_fractional_error) {
                    max_fractional_error = fractional_error;
                }
                if (fractional_error < min_fractional_error) {
                    min_fractional_error = fractional_error;
                }
                have_points = true;
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    if (!have_points) {
        max_fractional_error = 0.0;
        min_fractional_error = 0.0;
    }
#if DBGPRINT
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;
#endif //DBGPRINT

    // on the development machine, I get
    //    max_fractional_error = 0.0134406
    //    min_fractional_error = -0.00017682
    const double solution_tolerance = 2.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_scalar_field2, Spherical_bunch_fixture)
{
    // n.b. We don't shift frames here. We just want a beam that's spherical
    //      in the frame in which we are working.
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho)); // [C/m^3]
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper(), Commxx());

    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    bool have_points = false;
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain_sptr()->get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y + z * z);
                double phi_exact_ijk = gaussian_electric_potential(Q, r, sigma);
                phi_exact.get_grid_points()[i][j][k] = phi_exact_ijk;
                double phi_calc_ijk = phi2->get_grid_points()[i][j][k]
                        * phi2->get_normalization();
                double fractional_error = (phi_calc_ijk - phi_exact_ijk)
                        / phi_exact_ijk;
                if (fractional_error > max_fractional_error) {
                    max_fractional_error = fractional_error;
                }
                if (fractional_error < min_fractional_error) {
                    min_fractional_error = fractional_error;
                }
                have_points = true;
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    if (!have_points) {
        max_fractional_error = 0.0;
        min_fractional_error = 0.0;
    }
#if DBGPRINT
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;
#endif //DBGPRINT
    // on the development machine, I get (on one run)
    //        max_fractional_error = 0.0164322
    //        min_fractional_error = -0.0141608
    const double solution_tolerance = 3.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(extract_scalar_field, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho)); // [C/m^3]
    local_rho.reset();
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                BOOST_CHECK_CLOSE(
                        phi2->get_grid_points()[i][j][k]*phi2->get_normalization(),
                        phi->get_grid_points()[i][j][k]*phi->get_normalization(),
                        tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(extract_scalar_field_with_guards, Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper(), Commxx());
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    if (comm.get_rank() == 1) {
        std::cout << "phi->get_lower() = " << phi->get_lower() << std::endl;
        std::cout << "phi->get_upper() = " << phi->get_upper() << std::endl;
        std::cout << "phi->get_lower_guard() = " << phi->get_lower_guard()
                << std::endl;
        std::cout << "phi->get_upper_guard() = " << phi->get_upper_guard()
                << std::endl;
    }
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    double have_points = false;
    for (int i = phi->get_lower_guard(); i < phi->get_upper_guard(); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain_sptr()->get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y + z * z);
                double phi_exact_ijk = gaussian_electric_potential(Q, r, sigma);
                phi_exact.get_grid_points()[i][j][k] = phi_exact_ijk;
                double phi_calc_ijk = phi->get_grid_points()[i][j][k]
                        * phi->get_normalization();
                double fractional_error = (phi_calc_ijk - phi_exact_ijk)
                        / phi_exact_ijk;
                if (fractional_error > max_fractional_error) {
                    max_fractional_error = fractional_error;
                }
                if (fractional_error < min_fractional_error) {
                    min_fractional_error = fractional_error;
                }
                have_points = true;
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    if (!have_points) {
        max_fractional_error = 0.0;
        min_fractional_error = 0.0;
    }
#if DBGPRINT
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;
#endif

    // on the development machine, I get
    //    max_fractional_error = 0.0134406
    //    min_fractional_error = -0.00017682
    const double solution_tolerance = 2.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_local_electric_field_component_exact_rho,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    bool have_points = false;
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        double max_fractional_error = -2.0;
        double min_fractional_error = 2.0;
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    double z, y, x;
                    local_En->get_domain_sptr()->get_cell_coordinates(i, j, k,
                            z, y, x);
                    double r = std::sqrt(x * x + y * y + z * z);
                    double var;
                    if (component == 0) {
                        var = x;
                    } else if (component == 1) {
                        var = y;
                    } else if (component == 2) {
                        var = z;
                    }
                    double En_exact_ijk = gaussian_electric_field_component(Q,
                            r, sigma, var);
                    double En_calc_ijk = local_En->get_grid_points()[i][j][k]
                            * local_En->get_normalization();
                    double fractional_error = (En_calc_ijk - En_exact_ijk)
                            / En_exact_ijk;
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                    have_points = true;
                }
            }
        }
        if (!have_points) {
            max_fractional_error = 0.0;
            min_fractional_error = 0.0;
        }
#if DBGPRINT
	std::cout << "max_fractional_error = " << max_fractional_error
		  << std::endl;
	std::cout << "min_fractional_error = " << min_fractional_error
		  << std::endl;
#endif //DBGPRINT
        // on the development machine, I get
        //        max_fractional_error = 0.148878
        //        min_fractional_error = -0.0393951
        //        max_fractional_error = 0.0932745
        //        min_fractional_error = -0.0141639
        //        max_fractional_error = 0.067946
        //        min_fractional_error = -0.00695024

        const double field_tolerance[] = { 20.0e-2, 12.0e-2, 8.0e-2 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_electric_field_component_exact_rho_gatherv_bcast,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component_gatherv_bcast(
                        *local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_electric_field_component_exact_rho_allgatherv,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component_allgatherv(
                        *local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_electric_field_component_exact_rho_allreduce,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component_allreduce(
                        *local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_electric_field_component_exact_rho,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component(*local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_scalar_field2_G2linear_exact_rho, 
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper(), Commxx());
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    double have_points = false;
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain_sptr()->get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y + z * z);
                double phi_exact_ijk = gaussian_electric_potential(Q, r, sigma);
                phi_exact.get_grid_points()[i][j][k] = phi_exact_ijk;
                double phi_calc_ijk = phi2->get_grid_points()[i][j][k]
                        * phi2->get_normalization();
                double fractional_error = (phi_calc_ijk - phi_exact_ijk)
                        / phi_exact_ijk;
                if (fractional_error > max_fractional_error) {
                    max_fractional_error = fractional_error;
                }
                if (fractional_error < min_fractional_error) {
                    min_fractional_error = fractional_error;
                }
                have_points = true;
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    if (!have_points) {
        max_fractional_error = 0.0;
        min_fractional_error = 0.0;
    }
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;

    // on the development machine, I get
    //    max_fractional_error = 0.000546788
    //    min_fractional_error = -0.0216729
    const double solution_tolerance = 4.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_scalar_field2_G2linear_particles,
        Spherical_bunch_fixture)
{
    // n.b. We don't shift frames here. We just want a beam that's spherical
    //      in the frame in which we are working.
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho)); // [C/m^3]
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper(), Commxx());

    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    bool have_points = false;
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain_sptr()->get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y + z * z);
                double phi_exact_ijk = gaussian_electric_potential(Q, r, sigma);
                phi_exact.get_grid_points()[i][j][k] = phi_exact_ijk;
                double phi_calc_ijk = phi2->get_grid_points()[i][j][k]
                        * phi2->get_normalization();
                double fractional_error = (phi_calc_ijk - phi_exact_ijk)
                        / phi_exact_ijk;
                if (fractional_error > max_fractional_error) {
                    max_fractional_error = fractional_error;
                }
                if (fractional_error < min_fractional_error) {
                    min_fractional_error = fractional_error;
                }
                have_points = true;
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    if (!have_points) {
        max_fractional_error = 0.0;
        min_fractional_error = 0.0;
    }
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;

    // on the development machine, I get (on one run)
    //        max_fractional_error = 0.00863058
    //        min_fractional_error = -0.0375927
    const double solution_tolerance = 8.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(extract_scalar_field_G2linear,
        Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho)); // [C/m^3]
    local_rho.reset();
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                BOOST_CHECK_CLOSE(
                        phi2->get_grid_points()[i][j][k]*phi2->get_normalization(),
                        phi->get_grid_points()[i][j][k]*phi->get_normalization(),
                        tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(extract_scalar_field_G2linear_with_guards,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper(), Commxx());
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    if (comm.get_rank() == 1) {
        std::cout << "phi->get_lower() = " << phi->get_lower() << std::endl;
        std::cout << "phi->get_upper() = " << phi->get_upper() << std::endl;
        std::cout << "phi->get_lower_guard() = " << phi->get_lower_guard()
                << std::endl;
        std::cout << "phi->get_upper_guard() = " << phi->get_upper_guard()
                << std::endl;
    }
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    double have_points = false;
    for (int i = phi->get_lower_guard(); i < phi->get_upper_guard(); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain_sptr()->get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y + z * z);
                double phi_exact_ijk = gaussian_electric_potential(Q, r, sigma);
                phi_exact.get_grid_points()[i][j][k] = phi_exact_ijk;
                double phi_calc_ijk = phi->get_grid_points()[i][j][k]
                        * phi->get_normalization();
                double fractional_error = (phi_calc_ijk - phi_exact_ijk)
                        / phi_exact_ijk;
                if (fractional_error > max_fractional_error) {
                    max_fractional_error = fractional_error;
                }
                if (fractional_error < min_fractional_error) {
                    min_fractional_error = fractional_error;
                }
                have_points = true;
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    if (!have_points) {
        max_fractional_error = 0.0;
        min_fractional_error = 0.0;
    }
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;

    // on the development machine, I get
    //        max_fractional_error = 0.000546788
    //        min_fractional_error = -0.0216729
    const double solution_tolerance = 4.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_local_electric_field_component_G2linear_exact_rho,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    bool have_points = false;
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        double max_fractional_error = -2.0;
        double min_fractional_error = 2.0;
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    double z, y, x;
                    local_En->get_domain_sptr()->get_cell_coordinates(i, j, k,
                            z, y, x);
                    double r = std::sqrt(x * x + y * y + z * z);
                    double var;
                    if (component == 0) {
                        var = x;
                    } else if (component == 1) {
                        var = y;
                    } else if (component == 2) {
                        var = z;
                    }
                    double En_exact_ijk = gaussian_electric_field_component(Q,
                            r, sigma, var);
                    double En_calc_ijk = local_En->get_grid_points()[i][j][k]
                            * local_En->get_normalization();
                    double fractional_error = (En_calc_ijk - En_exact_ijk)
                            / En_exact_ijk;
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                    have_points = true;
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
        //        max_fractional_error = 0.144236
        //        min_fractional_error = -0.125599
        //        max_fractional_error = 0.0902098
        //        min_fractional_error = -0.0926057
        //        max_fractional_error = 0.0688483
        //        min_fractional_error = -0.0878926
        const double field_tolerance[] = { 2.0e-1, 2.0e-1, 2.0e-1 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}

BOOST_FIXTURE_TEST_CASE(
        get_global_electric_field_component_G2linear_exact_rho_gatherv_bcast,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component_gatherv_bcast(
                        *local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(
        get_global_electric_field_component_G2linear_exact_rho_allgatherv,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component_allgatherv(
                        *local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(
        get_global_electric_field_component_G2linear_exact_rho_allreduce,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component_allreduce(
                        *local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_electric_field_component_G2linear_exact_rho,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(comm, grid_shape);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_linear()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr En(
                space_charge.get_global_electric_field_component(*local_En)); // [V/m]
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain_sptr()->get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain_sptr()->get_grid_shape()[2]; ++k) {
                    BOOST_CHECK_CLOSE(local_En->get_grid_points()[i][j][k],
                            En->get_grid_points()[i][j][k], tolerance);
                }
            }
        }
        BOOST_CHECK_CLOSE(local_En->get_normalization(),
                En->get_normalization(), tolerance);
    }
}
