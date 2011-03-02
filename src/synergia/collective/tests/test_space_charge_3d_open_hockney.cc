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

BOOST_GLOBAL_FIXTURE(MPI_fixture)

const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.7e11;
const int total_num = 10000;
const double total_energy = 125.0;
struct Ellipsoidal_bunch_fixture
{
    Ellipsoidal_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm(MPI_COMM_WORLD), bunch(reference_particle,
                total_num, real_num, comm), distribution(0, comm),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup ellipsoidal bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        stdx = 1.1e-3;
        stdy = 2.3e-3;
        stdz = 3.5e-3;
        covariances[0][0] = stdx * stdx;
        covariances[2][2] = stdy * stdy;
        covariances[4][4] = stdz * stdz;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 16;
        grid_shape[1] = 24;
        grid_shape[2] = 32;
    }

    ~Ellipsoidal_bunch_fixture()
    {
        BOOST_TEST_MESSAGE("tear down ellipsoidal bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
    Random_distribution distribution;
    double stdx, stdy, stdz;
    std::vector<int > grid_shape;
};

struct Spherical_bunch_fixture
{
    Spherical_bunch_fixture() :
        four_momentum(mass, total_energy), reference_particle(charge,
                four_momentum), comm(MPI_COMM_WORLD), bunch(reference_particle,
                total_num, real_num, comm), distribution(0, comm),
                grid_shape(3)
    {
        BOOST_TEST_MESSAGE("setup Spherical bunch fixture");
        MArray2d covariances(boost::extents[6][6]);
        MArray1d means(boost::extents[6]);
        for (int i = 0; i < 6; ++i) {
            means[i] = 0.0;
            for (int j = i; j < 6; ++j) {
                covariances[i][j] = 0.0;
            }
        }
        sigma = 1.3e-3;
        covariances[0][0] = sigma * sigma;
        covariances[2][2] = sigma * sigma;
        covariances[4][4] = sigma * sigma;
        covariances[1][1] = covariances[3][3] = covariances[5][5] = 1.0;
        populate_6d(distribution, bunch, means, covariances);
        grid_shape[0] = 16;
        grid_shape[1] = 24;
        grid_shape[2] = 32;
    }

    ~Spherical_bunch_fixture()
    {
        BOOST_TEST_MESSAGE("tear down Spherical bunch fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
    Random_distribution distribution;
    double sigma;
    std::vector<int > grid_shape;
};

const double domain_min = -2.0;
const double domain_max = 2.0;
const double domain_offset = 0.0;
const int toy_grid_shape[] = { 4, 6, 8 };
const int toy_total_num = 1;
const double toy_real_num = 1.0e20;

struct Toy_bunch_fixture
{
    Toy_bunch_fixture() :
                four_momentum(mass, total_energy),
                reference_particle(pconstants::proton_charge, four_momentum),
                comm(MPI_COMM_WORLD),
                bunch(reference_particle, total_num, toy_real_num, comm),
                physical_size(3),
                physical_offset(3),
                cell_size(3),
                grid_shape(3),
                expected(
                        boost::extents[toy_grid_shape[2]][toy_grid_shape[1]][toy_grid_shape[0]])
    {
        for (int i = 0; i < 3; ++i) {
            grid_shape[i] = toy_grid_shape[i];
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            cell_size[i] = (domain_max - domain_min) / grid_shape[i];
        }
        for (unsigned int i = 0; i < expected.shape()[0]; ++i) {
            for (unsigned int j = 0; j < expected.shape()[1]; ++j) {
                for (unsigned int k = 0; k < expected.shape()[2]; ++k) {
                    expected[i][j][k] = 0.0;
                }
            }
        }
        density_norm = (toy_real_num / total_num) * pconstants::e
                / (cell_size[0] * cell_size[1] * cell_size[2]);
    }

    ~Toy_bunch_fixture()
    {
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
    double density_norm;

    std::vector<double > physical_size, physical_offset, cell_size;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray3d expected;
};

const double tolerance = 1.0e-12;

BOOST_AUTO_TEST_CASE(construct1)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx comm(MPI_COMM_WORLD);

    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
}

BOOST_AUTO_TEST_CASE(get_n_sigma)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 16;
    Commxx comm(MPI_COMM_WORLD);
    double my_n_sigma(7.0);

    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm, 0.0,
            my_n_sigma);
    BOOST_CHECK_CLOSE(space_charge.get_n_sigma(), my_n_sigma, tolerance);
}

BOOST_FIXTURE_TEST_CASE(update_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    space_charge.update_domain(bunch);
}

BOOST_FIXTURE_TEST_CASE(get_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    space_charge.update_domain(bunch);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[0],
            space_charge.get_n_sigma() * stdz, tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[1],
            space_charge.get_n_sigma() * stdy, tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[2],
            space_charge.get_n_sigma() * stdx, tolerance);
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[0],
                    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[1],
                    0.0, tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[2],
                    0.0, tolerance));
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[0],
            grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[1],
            grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[2],
            grid_shape[0]);
}

BOOST_FIXTURE_TEST_CASE(get_doubled_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    space_charge.update_domain(bunch);
    BOOST_CHECK_EQUAL(space_charge.get_doubled_domain_sptr()->get_grid_shape()[0],
            2*grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_doubled_domain_sptr()->get_grid_shape()[1],
            2*grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_doubled_domain_sptr()->get_grid_shape()[2],
            2*grid_shape[0]);
}

BOOST_FIXTURE_TEST_CASE(set_fixed_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    space_charge.update_domain(bunch);
    Rectangular_grid_domain_sptr domain_sptr = space_charge.get_domain_sptr();
    double scale = 1.23;
    double shift = 0.2345;
    std::vector<double > scaled_size(3), shifted_offset(3);
    for (int i = 0; i < 3; ++i) {
        scaled_size[i] = scale * domain_sptr->get_physical_size()[i];
        shifted_offset[i] = shift + domain_sptr->get_physical_offset()[i];
    }
    Rectangular_grid_domain_sptr new_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    scaled_size, shifted_offset, domain_sptr->get_grid_shape(),
                    domain_sptr->is_periodic()));
    space_charge.set_fixed_domain(new_domain_sptr);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[0],
            scaled_size[0], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[1],
            scaled_size[1], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[2],
            scaled_size[2], tolerance);
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[0],
                    shifted_offset[0], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[1],
                    shifted_offset[1], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[2],
                    shifted_offset[2], tolerance));
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[0],
            grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[1],
            grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[2],
            grid_shape[0]);

    // make sure that update domain has no effect
    space_charge.update_domain(bunch);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[0],
            scaled_size[0], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[1],
            scaled_size[1], tolerance);
    BOOST_CHECK_CLOSE(space_charge.get_domain_sptr()->get_physical_size()[2],
            scaled_size[2], tolerance);
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[0],
                    shifted_offset[0], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[1],
                    shifted_offset[1], tolerance));
    BOOST_CHECK(floating_point_equal(
                    space_charge.get_domain_sptr()->get_physical_offset()[2],
                    shifted_offset[2], tolerance));
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[0],
            grid_shape[2]);
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[1],
            grid_shape[1]);
    BOOST_CHECK_EQUAL(space_charge.get_domain_sptr()->get_grid_shape()[2],
            grid_shape[0]);
}

BOOST_FIXTURE_TEST_CASE(set_fixed_domain_bad_shape, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    space_charge.update_domain(bunch);
    Rectangular_grid_domain_sptr domain_sptr = space_charge.get_domain_sptr();
    std::vector<int > bad_shape(domain_sptr->get_grid_shape());
    // not bad yet
    Rectangular_grid_domain_sptr not_bad_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    space_charge.set_fixed_domain(not_bad_domain_sptr);

    bool caught_error;
    // bad in 0
    caught_error = false;
    ++bad_shape[0];
    Rectangular_grid_domain_sptr bad0_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    try {
        space_charge.set_fixed_domain(bad0_domain_sptr);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);

    // bad in 1
    caught_error = false;
    --bad_shape[0];
    ++bad_shape[1];
    Rectangular_grid_domain_sptr bad1_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    try {
        space_charge.set_fixed_domain(bad1_domain_sptr);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);

    // bad in 2
    caught_error = false;
    --bad_shape[1];
    ++bad_shape[2];
    Rectangular_grid_domain_sptr bad2_domain_sptr =
            Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    domain_sptr->get_physical_size(),
                    domain_sptr->get_physical_offset(), bad_shape,
                    domain_sptr->is_periodic()));
    try {
        space_charge.set_fixed_domain(bad2_domain_sptr);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(get_local_charge_density, Toy_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    std::vector<int > grid_shape_zyx(3);
    grid_shape_zyx[0] = grid_shape[2];
    grid_shape_zyx[1] = grid_shape[1];
    grid_shape_zyx[2] = grid_shape[0];
    Rectangular_grid_domain_sptr domain_sptr = Rectangular_grid_domain_sptr(
            new Rectangular_grid_domain(physical_size, physical_offset,
                    grid_shape_zyx, false));
    space_charge.set_fixed_domain(domain_sptr);
    std::vector<int > center_zyx(3);
    center_zyx[0] = grid_shape_zyx[0] / 2;
    center_zyx[1] = grid_shape_zyx[1] / 2;
    center_zyx[2] = grid_shape_zyx[2] / 2;

    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;
    Rectangular_grid_sptr local_charge_density(
            space_charge.get_local_charge_density(bunch));

    for (int i = center_zyx[0] - 1; i < center_zyx[0] + 1; ++i) {
        for (int j = center_zyx[1] - 1; j < center_zyx[1] + 1; ++j) {
            for (int k = center_zyx[2] - 1; k < center_zyx[2] + 1; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
            }
        }
    }
    multi_array_check_equal(local_charge_density->get_grid_points(), expected,
            100 * tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_global_charge_density2, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    Rectangular_grid_sptr local_rho = space_charge.get_local_charge_density(
            bunch); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2 =
            space_charge.get_global_charge_density2(*local_rho); // [C/m^3]
    std::vector<int > nondoubled_shape(
            local_rho->get_domain_sptr()->get_grid_shape());
    std::vector<int > doubled_shape(rho2->get_domain_sptr()->get_grid_shape());
    BOOST_CHECK_EQUAL(2*local_rho->get_domain_sptr()->get_grid_shape()[0],
            doubled_shape[0]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain_sptr()->get_grid_shape()[1],
            doubled_shape[1]);
    BOOST_CHECK_EQUAL(2*local_rho->get_domain_sptr()->get_grid_shape()[2],
            doubled_shape[2]);
    for (int i = 0; i < doubled_shape[0]; ++i) {
        for (int j = 0; j < doubled_shape[1]; ++j) {
            for (int k = 0; k < doubled_shape[2]; ++k) {
                if ((i < nondoubled_shape[0]) && (j < nondoubled_shape[1])
                        && (k < nondoubled_shape[2])) {
                    BOOST_CHECK_CLOSE(rho2->get_grid_points()[i][j][k],
                            local_rho->get_grid_points()[i][j][k], tolerance);
                } else {
                    BOOST_CHECK_EQUAL(rho2->get_grid_points()[i][j][k], 0.0);
                }
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_green_fn2_pointlike, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    space_charge.update_domain(bunch);
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike());
    MArray3d_ref G2_a(G2->get_grid_points());
    double norm = G2->get_normalization();
    int imirror, jmirror, kmirror;
    double dz, dy, dx;
    const double coeff = 2.8;
    double G000 = coeff / std::min(G2->get_domain_sptr()->get_cell_size()[0],
            std::min(G2->get_domain_sptr()->get_cell_size()[1],
                    G2->get_domain_sptr()->get_cell_size()[2]));

    int i_max = std::min(G2->get_upper(),
            G2->get_domain_sptr()->get_grid_shape()[0] / 2);
    for (int i = G2->get_lower(); i < i_max; ++i) {
        dz = i * G2->get_domain_sptr()->get_cell_size()[0];
        imirror = G2->get_domain_sptr()->get_grid_shape()[0] - i;
        if (imirror == G2->get_domain_sptr()->get_grid_shape()[0]) {
            imirror = i;
        }
        for (int j = 0; j < G2->get_domain_sptr()->get_grid_shape()[1] / 2; ++j) {
            dy = j * G2->get_domain_sptr()->get_cell_size()[1];
            jmirror = G2->get_domain_sptr()->get_grid_shape()[1] - j;
            if (jmirror == G2->get_domain_sptr()->get_grid_shape()[1]) {
                jmirror = j;
            }
            for (int k = 0; k < G2->get_domain_sptr()->get_grid_shape()[2] / 2; ++k) {
                dx = k * G2->get_domain_sptr()->get_cell_size()[2];
                kmirror = G2->get_domain_sptr()->get_grid_shape()[2] - k;
                if (kmirror == G2->get_domain_sptr()->get_grid_shape()[2]) {
                    kmirror = k;
                }
                double G;
                if ((i == 0) && (j == 0) && (k == 0)) {
                    G = G000;
                } else {
                    G = 1 / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
                BOOST_CHECK_CLOSE(G2_a[i][j][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[imirror][j][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[imirror][jmirror][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[imirror][jmirror][kmirror]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[imirror][j][kmirror]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][jmirror][k]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][jmirror][kmirror]*norm, G, tolerance);
                BOOST_CHECK_CLOSE(G2_a[i][j][kmirror]*norm, G, tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_green_fn2_no_domain, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    bool caught_error = false;
    try {
        Distributed_rectangular_grid_sptr G2(
                space_charge.get_green_fn2_pointlike());
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}

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
    for (int i = 0; i < doubled_shape[0]; ++i) {
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
    for (int i = 0; i < nondoubled_shape[0]; ++i) {
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
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper());

    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
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
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    //    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    //    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;

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
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    Rectangular_grid_sptr local_rho(space_charge.get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2(space_charge.get_global_charge_density2(
            *local_rho)); // [C/m^3]
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid phi_exact(phi2->get_domain_sptr(),
            phi2->get_lower(), phi2->get_upper());

    double Q = bunch.get_real_num() * bunch.get_particle_charge()
            * pconstants::e;
    std::vector<int > nondoubled_shape(
            space_charge.get_domain_sptr()->get_grid_shape());
    double max_fractional_error = -2.0;
    double min_fractional_error = 2.0;
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
                // BOOST_CHECK_CLOSE(phi_calc_ijk, phi_exact_ijk, solution_tolerance);
            }
        }
    }
    //    std::cout << "max_fractional_error = " << max_fractional_error << std::endl;
    //    std::cout << "min_fractional_error = " << min_fractional_error << std::endl;

    // on the development machine, I get (on one run)
    //        max_fractional_error = 0.0164322
    //        min_fractional_error = -0.0141608
    const double solution_tolerance = 3.0e-2;
    BOOST_CHECK(std::abs(max_fractional_error) < solution_tolerance);
    BOOST_CHECK(std::abs(min_fractional_error) < solution_tolerance);
}

BOOST_FIXTURE_TEST_CASE(extract_scalar_field, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
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

BOOST_FIXTURE_TEST_CASE(get_local_electric_field_component_exact_rho,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
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
    phi->fill_guards(comm);
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
                        var = z;
                    } else if (component == 1) {
                        var = y;
                    } else if (component == 2) {
                        var = x;
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
                }
            }
        }
        //std::cout << "max_fractional_error = " << max_fractional_error
        //        << std::endl;
        //std::cout << "min_fractional_error = " << min_fractional_error
        //        << std::endl;
        // on the development machine, I get
        //    max_fractional_error = 0.067946
        //    min_fractional_error = -0.00695024
        //    max_fractional_error = 0.0932745
        //    min_fractional_error = -0.0141639
        //    max_fractional_error = 0.148878
        //    min_fractional_error = -0.0393951
        const double field_tolerance[] = { 8.0e-2, 12.0e-2, 20.0e-2 };
        BOOST_CHECK(std::abs(max_fractional_error) < field_tolerance[component]);
        BOOST_CHECK(std::abs(min_fractional_error) < field_tolerance[component]);
    }
}

BOOST_FIXTURE_TEST_CASE(get_global_electric_field_component_exact_rho,
        Spherical_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    Distributed_rectangular_grid_sptr rho2(get_gaussian_rho2(space_charge,
            bunch, sigma));
    Distributed_rectangular_grid_sptr
            G2(space_charge.get_green_fn2_pointlike()); // [1/m]
    Distributed_rectangular_grid_sptr phi2(space_charge.get_scalar_field2(
            *rho2, *G2)); // [V]
    Distributed_rectangular_grid_sptr phi(space_charge.extract_scalar_field(
            *phi2));
    phi->fill_guards(comm);
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
    }
}

BOOST_FIXTURE_TEST_CASE(apply, Ellipsoidal_bunch_fixture)
{
    Space_charge_3d_open_hockney space_charge(grid_shape, false, comm);
    Step dummy_step(1.0);
    space_charge.apply(bunch, 1.0, dummy_step);
    // jfa : n.b. this test is incomplete
}
