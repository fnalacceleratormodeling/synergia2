#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/bunch.h"
#include "components/collective/deposit.h"
#include "components/foundation/physical_constants.h"
#include "components/bunch/bunch.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

void
compare_multi_array(Const_MArray3d_ref const& a, Const_MArray3d_ref const& b,
        double tolerance)
{
    BOOST_CHECK_EQUAL(a.shape()[0],b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1],b.shape()[1]);
    BOOST_CHECK_EQUAL(a.shape()[2],b.shape()[2]);
    for (int i = 0; i < a.shape()[0]; ++i) {
        for (int j = 0; j < a.shape()[1]; ++j) {
            for (int k = 0; k < a.shape()[2]; ++k) {
                if (a[i][j][k] == 0.0) {
                    BOOST_CHECK_SMALL(b[i][j][k], tolerance);
                } else {
                    BOOST_CHECK_CLOSE(a[i][j][k], b[i][j][k], tolerance);
                }
            }
        }
    }
}

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 10;
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
                constants::proton_charge, four_momentum), comm(MPI_COMM_WORLD),
                bunch(reference_particle, total_num, real_num, comm),
                physical_size(3), physical_offset(3), grid_shape(3), expected(
                        boost::extents[grid_size][grid_size][grid_size])
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            grid_shape[i] = grid_size;
        }
        for (int i = 0; i < expected.shape()[0]; ++i) {
            for (int j = 0; j < expected.shape()[1]; ++j) {
                for (int k = 0; k < expected.shape()[2]; ++k) {
                    expected[i][j][k] = 0.0;
                }
            }
        }
    }

    ~Fixture()
    {
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Rectangular_grid_sptr rho_grid_sptr;
    MArray3d expected;
};

BOOST_FIXTURE_TEST_CASE(no_particles, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(0);
    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(origin_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 0.125;
            }
        }
    }

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(xedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = domain_min;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    expected[1][1][0] = 0.125;
    expected[1][2][0] = 0.125;
    expected[2][1][0] = 0.125;
    expected[2][2][0] = 0.125;

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(yedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = domain_min;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    expected[1][0][1] = 0.125;
    expected[1][0][2] = 0.125;
    expected[2][0][1] = 0.125;
    expected[2][0][2] = 0.125;

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(zedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_min;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.125;
    expected[0][1][2] = 0.125;
    expected[0][2][1] = 0.125;
    expected[0][2][2] = 0.125;

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(xrightedge_particle, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = domain_max;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    expected[1][1][3] = 0.125;
    expected[1][2][3] = 0.125;
    expected[2][1][3] = 0.125;
    expected[2][2][3] = 0.125;

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}


BOOST_FIXTURE_TEST_CASE(zedge_particle_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_min;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.125;
    expected[0][1][2] = 0.125;
    expected[0][2][1] = 0.125;
    expected[0][2][2] = 0.125;

    expected[3][1][1] = 0.125;
    expected[3][1][2] = 0.125;
    expected[3][2][1] = 0.125;
    expected[3][2][2] = 0.125;

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(zrightedge_particle_periodic, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, true));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = domain_max;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);

    expected[0][1][1] = 0.125;
    expected[0][1][2] = 0.125;
    expected[0][2][1] = 0.125;
    expected[0][2][2] = 0.125;

    expected[3][1][1] = 0.125;
    expected[3][1][2] = 0.125;
    expected[3][2][1] = 0.125;
    expected[3][2][2] = 0.125;

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}

BOOST_FIXTURE_TEST_CASE(additive_deposit, Fixture)
{
    rho_grid_sptr = Rectangular_grid_sptr(new Rectangular_grid(physical_size,
            physical_offset, grid_shape, false));
    bunch.set_local_num(1);
    bunch.get_local_particles()[0][0] = 0;
    bunch.get_local_particles()[0][2] = 0;
    bunch.get_local_particles()[0][4] = 0;

    deposit_charge_rectangular(*rho_grid_sptr, bunch);
    deposit_charge_rectangular(*rho_grid_sptr, bunch, false);

    for (int i = 1; i < 3; ++i) {
        for (int j = 1; j < 3; ++j) {
            for (int k = 1; k < 3; ++k) {
                expected[i][j][k] = 2*0.125;
            }
        }
    }

    compare_multi_array(rho_grid_sptr->get_grid_points(), expected, tolerance);
}
