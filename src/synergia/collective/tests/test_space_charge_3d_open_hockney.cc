#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/floating_point.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const int charge = pconstants::proton_charge;
const double mass = pconstants::mp;
const double real_num = 1.7e11;
const int total_num = 1000;
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

