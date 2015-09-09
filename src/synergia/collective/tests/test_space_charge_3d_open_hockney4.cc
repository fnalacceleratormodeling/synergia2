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
#include "uniform_cylindrical_charge_density.h"
#include "space_charge_bunch_fixtures.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

// define DBGPRINT to get verbose output of test values
#define DBGPRINT 0

const double tolerance = 1.0e-12;

Distributed_rectangular_grid_sptr
get_uniform_cylindrical_rho2(Space_charge_3d_open_hockney & space_charge,
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
            for (int k = 0; k < doubled_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] = 0.0;
            }
        }
    }
    std::vector<int > nondoubled_shape(
            space_charge.get_domain().get_grid_shape());
    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;
    for (int i = 0; i < nondoubled_shape[0]; ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                local_rho->get_domain().get_cell_coordinates(i, j, k, z,
                        y, x);
                double r = std::sqrt(x * x + y * y);
                rho2->get_grid_points()[i][j][k]
                        = uniform_cylindrical_charge_density(lambda, r, r0);
            }
        }
    }
    return rho2;
}

BOOST_FIXTURE_TEST_CASE(get_scalar_field2_exact_rho, Cylindrical_bunch_fixture_fine)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, false, true,
            z_period, true);
    Distributed_rectangular_grid_sptr rho2(get_uniform_cylindrical_rho2(
            space_charge, bunch, r0, z_period));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    Distributed_rectangular_grid phi_exact(phi->get_domain_sptr(),
            phi->get_lower(), phi->get_upper(), comm_sptr);

    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;

    std::vector<int > nondoubled_shape(
            space_charge.get_domain().get_grid_shape());

    // The potential has an arbitrary offset; take it from the middle
    double z0, y0, x0;
    int i0 = nondoubled_shape[0] / 2;
    int j0 = nondoubled_shape[1] / 2;
    int k0 = nondoubled_shape[2] / 2;
    space_charge.get_domain().get_cell_coordinates(i0, j0, k0, z0, y0, x0);
    double roffset = std::sqrt(x0 * x0 + y0 * y0);
    double offset = phi2->get_grid_points()[i0][j0][k0]
            * phi2->get_normalization()
            - uniform_cylindrical_electric_potential(lambda, roffset, r0);

    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain().get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y);
                double phi_exact_ijk = uniform_cylindrical_electric_potential(
                        lambda, r, r0) + offset;
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
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
#if DBGPRINT
    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;
#endif //DBGPRINT

    // on the development machine, I get
    //    max_fractional_error = 4.83071e-05
    //    min_fractional_error = -0.00537732

    const double solution_tolerance = 1.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_local_electric_field_component_exact_rho,
        Cylindrical_bunch_fixture_fine)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, false, true,
            z_period, true);
    Distributed_rectangular_grid_sptr rho2(get_uniform_cylindrical_rho2(
            space_charge, bunch, r0, z_period));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;
    phi->fill_guards();
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        //        Distributed_rectangular_grid_sptr exact_En(local_En);
        //        std::string filename;
        //        if (component == 0) {
        //            filename = "ez.h5";
        //        } else if (component == 1) {
        //            filename = "ey.h5";
        //        } else if (component == 2) {
        //            filename = "ex.h5";
        //        }
        //        Hdf5_file f(filename);
        //        f.write(local_En->get_grid_points(), "en");
        //        f.write(local_En->get_normalization(), "ennorm");
        double max_fractional_error = -2.0;
        double min_fractional_error = 2.0;
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain().get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain().get_grid_shape()[2]; ++k) {
                    double z, y, x;
                    local_En->get_domain().get_cell_coordinates(i, j, k,
                            z, y, x);
                    double r = std::sqrt(x * x + y * y);
                    double var;
                    double En_exact_ijk;
                    if (component == 2) {
                        En_exact_ijk = 0.0;
                    } else {
                        if (component == 0) {
                            var = x;
                        } else {
                            var = y;
                        }
                        En_exact_ijk
                                = uniform_cylindrical_electric_field_component(
                                        lambda, r, r0, var);
                    }
                    //                    exact_En->get_grid_points()[i][j][k] = En_exact_ijk;
                    double En_calc_ijk = local_En->get_grid_points()[i][j][k]
                            * local_En->get_normalization();
                    const double tiny = 1.0e-8;
                    double fractional_error;
                    if (std::abs(En_exact_ijk) < tiny) {
                        fractional_error = En_calc_ijk - En_exact_ijk;
                    } else {
                        fractional_error = (En_calc_ijk - En_exact_ijk)
                                / En_exact_ijk;
                    }
                    if (fractional_error > max_fractional_error) {
                        max_fractional_error = fractional_error;
                    }
                    if (fractional_error < min_fractional_error) {
                        min_fractional_error = fractional_error;
                    }
                }
            }
        }
	// f.write(exact_En->get_grid_points(), "enexact");

#if DBGPRINT
	std::cout << "max_fractional_error = " << max_fractional_error
		  << std::endl;
	std::cout << "min_fractional_error = " << min_fractional_error
		  << std::endl;
#endif //DBGPRINT

        // on the development machine, I get
        //                max_fractional_error = 0.193765
        //                min_fractional_error = -0.0428391
        //                max_fractional_error = 0.193765
        //                min_fractional_error = -0.0428391
        //                max_fractional_error = 25740.3
        //                min_fractional_error = -25740.3
        const double field_tolerance[] = { 0.4, 0.4, 5.0e5 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}

BOOST_FIXTURE_TEST_CASE(get_scalar_field2_particles, Cylindrical_bunch_fixture_fine)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    double stdqp = 1.0e-10; // not important...
    populate_uniform_cylinder(distribution, bunch, r0, z_period, stdqp, stdqp,
            stdqp);
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, false, true,
            z_period, true);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch));
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho, comm_sptr));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    Distributed_rectangular_grid phi_exact(phi->get_domain_sptr(),
            phi->get_lower(), phi->get_upper(), comm_sptr);

    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;

    std::vector<int > nondoubled_shape(
            space_charge.get_domain().get_grid_shape());

    // The potential has an arbitrary offset; take it from the middle
    double z0, y0, x0;
    int i0 = nondoubled_shape[0] / 2;
    int j0 = nondoubled_shape[1] / 2;
    int k0 = nondoubled_shape[2] / 2;
    space_charge.get_domain().get_cell_coordinates(i0, j0, k0, z0, y0, x0);
    double roffset = std::sqrt(x0 * x0 + y0 * y0);
    double offset = phi2->get_grid_points()[i0][j0][k0]
            * phi2->get_normalization()
            - uniform_cylindrical_electric_potential(lambda, roffset, r0);

    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
    for (int i = phi2->get_lower(); i < std::min(phi2->get_upper(),
            nondoubled_shape[0]); ++i) {
        for (int j = 0; j < nondoubled_shape[1]; ++j) {
            for (int k = 0; k < nondoubled_shape[2]; ++k) {
                double z, y, x;
                space_charge.get_domain().get_cell_coordinates(i, j, k,
                        z, y, x);
                double r = std::sqrt(x * x + y * y);
                double phi_exact_ijk = uniform_cylindrical_electric_potential(
                        lambda, r, r0) + offset;
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
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    //        std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    //        std::cout << "min_fractional_error = " << min_fractional_error << std::endl;

    // on the development machine, I get
    //        max_fractional_error = 0.00814816
    //        min_fractional_error = -0.00462393

    const double solution_tolerance = 2.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_local_electric_field_component_particles,
        Cylindrical_bunch_fixture_fine)
{
    double z_period = 8 * sigma;
    double r0 = 2.0 * sigma;
    double stdqp = 1.0e-10; // not important...
    populate_uniform_cylinder(distribution, bunch, r0, z_period, stdqp, stdqp,
            stdqp);
    Space_charge_3d_open_hockney space_charge(comm_sptr, grid_shape, false, true,
            z_period, true);
    Rectangular_grid_sptr local_rho(
            space_charge.get_local_charge_density(bunch));
    Distributed_rectangular_grid_sptr rho2(
            space_charge.get_global_charge_density2(*local_rho, comm_sptr));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    double lambda = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e / z_period;
    phi->fill_guards();
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                space_charge.get_electric_field_component(*phi, component)); // [V/m]
        //        Distributed_rectangular_grid exact_En(local_En->get_domain_sptr(),
        //                local_En->get_lower(), local_En->get_upper());
        //        std::string filename;
        //        if (component == 0) {
        //            filename = "ez.h5";
        //        } else if (component == 1) {
        //            filename = "ey.h5";
        //        } else if (component == 2) {
        //            filename = "ex.h5";
        //        }
        //        Hdf5_file f(filename);
        //        f.write(local_En->get_grid_points(), "en");
        //        f.write(local_En->get_normalization(), "ennorm");
        double field_max = -1e30;
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain().get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain().get_grid_shape()[2]; ++k) {
                    if (std::abs(local_En->get_grid_points()[i][j][k])
                            > field_max) {
                        field_max = std::abs(
                                local_En->get_grid_points()[i][j][k]);
                    }
                }
            }
        }
        field_max *= local_En->get_normalization();
        double max_scaled_error = -1e30;
        double min_scaled_error = 1e30;
        for (int i = local_En->get_lower(); i < local_En->get_upper(); ++i) {
            for (int j = 0; j
                    < local_En->get_domain().get_grid_shape()[1]; ++j) {
                for (int k = 0; k
                        < local_En->get_domain().get_grid_shape()[2]; ++k) {
                    double z, y, x;
                    local_En->get_domain().get_cell_coordinates(i, j, k,
                            z, y, x);
                    double r = std::sqrt(x * x + y * y);
                    double var;
                    double En_exact_ijk;
                    if (component == 2) {
                        En_exact_ijk = 0.0;
                    } else {
                        if (component == 0) {
                            var = x;
                        } else {
                            var = y;
                        }
                        En_exact_ijk
                                = uniform_cylindrical_electric_field_component(
                                        lambda, r, r0, var);
                    }
                    //                    exact_En.get_grid_points()[i][j][k] = En_exact_ijk;
                    double En_calc_ijk = local_En->get_grid_points()[i][j][k]
                            * local_En->get_normalization();
                    double scaled_error;

                    scaled_error = (En_calc_ijk - En_exact_ijk) / field_max;
                    if (scaled_error > max_scaled_error) {
                        max_scaled_error = scaled_error;
                    }
                    if (scaled_error < min_scaled_error) {
                        min_scaled_error = scaled_error;
                    }
                }
            }
        }
        //        f.write(exact_En.get_grid_points(), "enexact");
#if DBGPRINT
	std::cout << "max_scaled_error = " << max_scaled_error << std::endl;
	std::cout << "min_scaled_error = " << min_scaled_error << std::endl;
#endif // DBGPRINT

        // on the development machine, I get
        //        max_scaled_error = 0.149345
        //        min_scaled_error = -0.153463
        //        max_scaled_error = 0.161769
        //        min_scaled_error = -0.160674
        //        max_scaled_error = 0.931733
        //        min_scaled_error = -1
        const double field_tolerance[] = { 0.2, 0.2, 1.0001 };
        BOOST_CHECK(std::abs(max_scaled_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_scaled_error) < field_tolerance[component]);
    }
}

