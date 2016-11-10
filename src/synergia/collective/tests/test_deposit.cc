#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/bunch.h"
#include "synergia/collective/deposit.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 21;
const double real_num = 2.0e12;
const double tolerance = 1.0e-12;
const double domain_min = -2.0;
const double domain_max = 2.0;
const double domain_offset = 0.0;
const int grid_size = 4;

struct Fixture
{
    Fixture() :
        four_momentum(mass, total_energy), reference_particle(
                pconstants::proton_charge, four_momentum),
                comm_sptr(new Commxx), bunch(reference_particle, total_num,
                        real_num, comm_sptr), physical_size(3), physical_offset(3),
                grid_shape(3), expected(
                        boost::extents[grid_size][grid_size][grid_size]),
                expected_2dc(boost::extents[grid_size][grid_size]),
                expected_1d(boost::extents[grid_size])
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            grid_shape[i] = grid_size;
        }
        for (unsigned int i = 0; i < expected.shape()[0]; ++i) {
            for (unsigned int j = 0; j < expected.shape()[1]; ++j) {
                for (unsigned int k = 0; k < expected.shape()[2]; ++k) {
                    expected[i][j][k] = 0.0;
                }
                expected_2dc[i][j] = 0.0;
            }
            expected_1d[i] = 0.0;
        }
        double cell_size = (domain_max - domain_min) / grid_size;
        density_norm = (real_num / total_num) * pconstants::e / (cell_size
                * cell_size * cell_size);
    }

    ~Fixture()
    {
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;
    double density_norm;

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray3d expected;
    MArray2dc expected_2dc;
    MArray1d  expected_1d;
};

BOOST_FIXTURE_TEST_CASE(no_particles, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(0);
    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(origin_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(in_domain, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    double hy= rho_grid_sptr->get_domain().get_cell_size()[1];
    bunch.set_local_num(2);
    bunch.get_local_particles()[0][0]=0.;
    bunch.get_local_particles()[0][2] =  physical_offset[1]+(grid_shape[1]/2*hy-0.49991*hy); //particle close to the left physicsl edge

  //  rho_grid_sptr->get_domain().get_left()[1]+physical_size[1]-0.6*hy;
   //
    bunch.get_local_particles()[0][4] = 0.;
    bunch.get_local_particles()[1][0]=0.;
    bunch.get_local_particles()[1][2] = physical_offset[1]-(grid_shape[1]/2*hy-0.49991*hy);//particle close to the right physicsl edge //rho_grid_sptr->get_domain().get_left()[1]+0.6*hy;

    bunch.get_local_particles()[1][4] = 0.;
    int ix, iy, iz;
    double offx, offy, offz;
    bool in_domain_r, in_domain_l;
    Const_MArray2d_ref parts(bunch.get_local_particles());
    in_domain_r=rho_grid_sptr->get_domain().get_leftmost_indices_offsets(
                    parts[0][4], parts[0][2], parts[0][0], iz, iy, ix, offz,
                    offy, offx);
  //  std::cout<<" right index iy="<<iy<<std::endl;
    in_domain_l=rho_grid_sptr->get_domain().get_leftmost_indices_offsets(
                    parts[1][4], parts[1][2], parts[1][0], iz, iy, ix, offz,
                    offy, offx);
   // std::cout<<" left index iy="<<iy<<std::endl;
    BOOST_CHECK (!in_domain_r);
    BOOST_CHECK (!in_domain_l);



}


BOOST_FIXTURE_TEST_CASE(x_displaced_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0]
            = rho_grid_sptr->get_domain().get_cell_size()[2];
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 2; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
            }
            expected[i][j][3] = 0.; //deposit on the grid's edge should be zero
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(x_displaced_particle1, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0]
            = -rho_grid_sptr->get_domain().get_cell_size()[2];
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 2; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
            }
            expected[i][j][0] = 0.; //deposit on the grid's edge should be zero
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(y_displaced_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2]
            = rho_grid_sptr->get_domain().get_cell_size()[1];
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 2; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
                expected[i][3][k] = 0.;//deposit on the grid's edge should be zero
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(y_displaced_particle1, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2]
            = -rho_grid_sptr->get_domain().get_cell_size()[1];
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 2; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
                expected[i][0][k] = 0.;//deposit on the grid's edge should be zero
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}


BOOST_FIXTURE_TEST_CASE(z_displaced_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4]
            = rho_grid_sptr->get_domain().get_cell_size()[0];

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 2; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
                expected[3][j][k] = 0.;//deposit on the grid's edge should be zer
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(z_displaced_particle1, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4]
            = -rho_grid_sptr->get_domain().get_cell_size()[0];

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 1; i < 2; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm;
                expected[0][j][k] = 0.;//deposit on the grid's edge should be zer
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}


BOOST_FIXTURE_TEST_CASE(z_displaced_particles_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    const int left_periods = 10;
    const int periods = 2 * left_periods + 1;
    bunch.set_local_num(periods);
    int n = 0;
    for (int period = -left_periods; period < left_periods + 1; ++period) {
        bunch.get_local_particles()[n][0] = 0;
        bunch.get_local_particles()[n][2] = 0;
        bunch.get_local_particles()[n][4]
                = rho_grid_sptr->get_domain().get_cell_size()[0] + period
                        * physical_size[2];
        ++n;
    }

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    for (int i = 2; i < 4; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125 * density_norm * periods;
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(xedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = domain_min;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[1][1][0] = 0.;//0.125 * density_norm;
    expected[1][2][0] = 0.;//0.125 * density_norm;
    expected[2][1][0] = 0.;//0.125 * density_norm;
    expected[2][2][0] = 0.;//0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(xedge_particle_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = domain_min;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[1][1][0] = 0.;//0.125 * density_norm;
    expected[1][2][0] = 0.;//0.125 * density_norm;
    expected[2][1][0] = 0.;//0.125 * density_norm;
    expected[2][2][0] = 0.;//0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(yedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = domain_min;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[1][0][1] = 0.;//0.125 * density_norm;
    expected[1][0][2] = 0.;//0.125 * density_norm;
    expected[2][0][1] = 0.;//0.125 * density_norm;
    expected[2][0][2] = 0.;//0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(yedge_particle_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = domain_min;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[1][0][1] = 0.;//0.125 * density_norm;
    expected[1][0][2] = 0.;//0.125 * density_norm;
    expected[2][0][1] = 0.;//0.125 * density_norm;
    expected[2][0][2] = 0.;//0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}
BOOST_FIXTURE_TEST_CASE(zedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_min;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.;// 0.125 * density_norm;
    expected[0][1][2] = 0.;//0.125 * density_norm;
    expected[0][2][1] = 0.;//0.125 * density_norm;
    expected[0][2][2] = 0.;//0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}



BOOST_FIXTURE_TEST_CASE(xrightedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = domain_max;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[1][1][3] = 0.;//0.125 * density_norm;
    expected[1][2][3] = 0.;//0.125 * density_norm;
    expected[2][1][3] = 0.;//0.125 * density_norm;
    expected[2][2][3] = 0.;//0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zedge_particle_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_min;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.125 * density_norm;
    expected[0][1][2] = 0.125 * density_norm;
    expected[0][2][1] = 0.125 * density_norm;
    expected[0][2][2] = 0.125 * density_norm;

    expected[3][1][1] = 0.125 * density_norm;
    expected[3][1][2] = 0.125 * density_norm;
    expected[3][2][1] = 0.125 * density_norm;
    expected[3][2][2] = 0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zedge_particles_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    const int left_periods = 10;
    const int periods = 2 * left_periods + 1;
    bunch.set_local_num(periods);
    int n = 0;
    for (int period = -left_periods; period < left_periods + 1; ++period) {
        bunch.get_local_particles()[n][0] = 0;
        bunch.get_local_particles()[n][2] = 0;
        bunch.get_local_particles()[n][4] = domain_min + period
                * physical_size[2];
        ++n;
    }

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.125 * density_norm * periods;
    expected[0][1][2] = 0.125 * density_norm * periods;
    expected[0][2][1] = 0.125 * density_norm * periods;
    expected[0][2][2] = 0.125 * density_norm * periods;

    expected[3][1][1] = 0.125 * density_norm * periods;
    expected[3][1][2] = 0.125 * density_norm * periods;
    expected[3][2][1] = 0.125 * density_norm * periods;
    expected[3][2][2] = 0.125 * density_norm * periods;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zrightedge_particle_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_max;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.125 * density_norm;
    expected[0][1][2] = 0.125 * density_norm;
    expected[0][2][1] = 0.125 * density_norm;
    expected[0][2][2] = 0.125 * density_norm;

    expected[3][1][1] = 0.125 * density_norm;
    expected[3][1][2] = 0.125 * density_norm;
    expected[3][2][1] = 0.125 * density_norm;
    expected[3][2][2] = 0.125 * density_norm;

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(additive_deposit, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch);
    deposit_charge_rectangular_zyx(*rho_grid_sptr, bunch, false);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 2 * 0.125 * density_norm;
            }
        }
    }

    multi_array_check_equal(rho_grid_sptr->get_grid_points(), expected,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(origin_particle_2d, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular_2d(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            expected_2dc[i][j] = 0.25 * density_norm;
        }
    }

    expected_1d[1] = 0.5;
    expected_1d[2] = 0.5;

    multi_complex_array_check_equal(rho_grid_sptr->get_grid_points_2dc(), expected_2dc,
            tolerance);

    multi_array_check_equal(rho_grid_sptr->get_grid_points_1d(), expected_1d,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zedge_particle_2d, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_min;

    deposit_charge_rectangular_2d(*rho_grid_sptr, bunch);

    expected_2dc[1][1] = 0.25 * density_norm;
    expected_2dc[1][2] = 0.25 * density_norm;
    expected_2dc[2][1] = 0.25 * density_norm;
    expected_2dc[2][2] = 0.25 * density_norm;

    expected_1d[0] = 0.5;

    multi_complex_array_check_equal(rho_grid_sptr->get_grid_points_2dc(), expected_2dc,
            tolerance);

    multi_array_check_equal(rho_grid_sptr->get_grid_points_1d(), expected_1d,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zrightedge_particle_2d, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_max;

    deposit_charge_rectangular_2d(*rho_grid_sptr, bunch);

    expected_2dc[1][1] = 0.25 * density_norm;
    expected_2dc[1][2] = 0.25 * density_norm;
    expected_2dc[2][1] = 0.25 * density_norm;
    expected_2dc[2][2] = 0.25 * density_norm;

    expected_1d[3] = 0.5;

    multi_complex_array_check_equal(rho_grid_sptr->get_grid_points_2dc(), expected_2dc,
            tolerance);

    multi_array_check_equal(rho_grid_sptr->get_grid_points_1d(), expected_1d,
            tolerance);
}


BOOST_FIXTURE_TEST_CASE(origin_particle_2d_bin, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    Raw_MArray2d bin(boost::extents[bunch.get_local_num()][6]);
    deposit_charge_rectangular_2d(*rho_grid_sptr, bin, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            expected_2dc[i][j] = 0.25 * density_norm;
        }
    }

    expected_1d[1] = 0.5;
    expected_1d[2] = 0.5;

    multi_complex_array_check_equal(rho_grid_sptr->get_grid_points_2dc(), expected_2dc,
            tolerance);

    multi_array_check_equal(rho_grid_sptr->get_grid_points_1d(), expected_1d,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zedge_particle_2d_bin, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_min;

    Raw_MArray2d bin(boost::extents[bunch.get_local_num()][6]);
    deposit_charge_rectangular_2d(*rho_grid_sptr, bin, bunch);

    expected_2dc[1][1] = 0.25 * density_norm;
    expected_2dc[1][2] = 0.25 * density_norm;
    expected_2dc[2][1] = 0.25 * density_norm;
    expected_2dc[2][2] = 0.25 * density_norm;

    expected_1d[0] = 0.5;

    multi_complex_array_check_equal(rho_grid_sptr->get_grid_points_2dc(), expected_2dc,
            tolerance);

    multi_array_check_equal(rho_grid_sptr->get_grid_points_1d(), expected_1d,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(zrightedge_particle_2d_bin, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_max;

    Raw_MArray2d bin(boost::extents[bunch.get_local_num()][6]);
    deposit_charge_rectangular_2d(*rho_grid_sptr, bin, bunch);

    expected_2dc[1][1] = 0.25 * density_norm;
    expected_2dc[1][2] = 0.25 * density_norm;
    expected_2dc[2][1] = 0.25 * density_norm;
    expected_2dc[2][2] = 0.25 * density_norm;

    expected_1d[3] = 0.5;

    multi_complex_array_check_equal(rho_grid_sptr->get_grid_points_2dc(), expected_2dc,
            tolerance);

    multi_array_check_equal(rho_grid_sptr->get_grid_points_1d(), expected_1d,
            tolerance);
}

