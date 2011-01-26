#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/space_charge_2d_bassetti_erskine.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/floating_point.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/hdf5_writer.h"
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

BOOST_AUTO_TEST_CASE(construct)
{
    Space_charge_2d_bassetti_erskine space_charge;
}

BOOST_FIXTURE_TEST_CASE(apply, Ellipsoidal_bunch_fixture)
{
    Space_charge_2d_bassetti_erskine space_charge;
    Step dummy_step;
    space_charge.apply(bunch, 1.0, dummy_step);
}
