#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/collective/impedance.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/step.h"
#include "synergia/simulation/operator.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-11;


double const  orbit_length=160.;
double const bunch_spacing=4.5;
int const  zgrid=40;
std::string const pipe_symmetry("x_parallel");


const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 100;
const double real_num = 2.0e12;
const double step_length = 1.23;
struct Fixture
{
    Fixture() :
        four_momentum(mass, total_energy), reference_particle(
                pconstants::proton_charge, four_momentum),
                comm(MPI_COMM_WORLD),
                 bunch(reference_particle, total_num,
                        real_num, comm), step(step_length)
    {
        BOOST_TEST_MESSAGE("setup fixture");
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
    Step step;
};

void
dummy_populate(Bunch &bunch, int offset = 0)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // coordinates
        for (int i = 0; i < 6; i += 2) {
            bunch.get_local_particles()[part][i] = 10.0 * (part + offset) + i;
        }
        // momenta
        for (int i = 1; i < 6; i += 2) {
            bunch.get_local_particles()[part][i] = 1e-4 * (10.0 * (part
                    + offset) + i);
        }
    }
}




BOOST_AUTO_TEST_CASE(test_constructor)
{
Impedance imped("test_wake.dat", orbit_length, bunch_spacing,zgrid, "circular",10);
}

BOOST_AUTO_TEST_CASE(wake_reading_circular)
{
Impedance imped("test_wake_cc.dat", orbit_length, bunch_spacing, zgrid, "circular",10); // three columns file
BOOST_CHECK_CLOSE(imped.get_z_coord()[0], 0.00614575, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_coord()[1], 0.0136406, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_coord()[2], 0.0241333, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_coord()[3], 0.037624, tolerance);

BOOST_CHECK_CLOSE(imped.get_z_wake()[0],  9.30907e+08, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_wake()[1],  5.79371e+08, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_wake()[2],  3.78051e+08 , tolerance);
BOOST_CHECK_CLOSE(imped.get_z_wake()[3],  2.55723e+08, tolerance);

BOOST_CHECK_CLOSE(imped.get_x_wake()[0],   -7.40775e+10   , tolerance);
BOOST_CHECK_CLOSE(imped.get_x_wake()[1],   -8.2791e+10, tolerance);
BOOST_CHECK_CLOSE(imped.get_x_wake()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(imped.get_x_wake()[3],   -9.39278e+10  , tolerance);

BOOST_CHECK_CLOSE(imped.get_y_wake()[0],   -7.40775e+10   , tolerance);
BOOST_CHECK_CLOSE(imped.get_y_wake()[1],   -8.2791e+10, tolerance);
BOOST_CHECK_CLOSE(imped.get_y_wake()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(imped.get_y_wake()[3],   -9.39278e+10  , tolerance);

}

BOOST_AUTO_TEST_CASE(wake_reading_parallel)
{
Impedance imped("test_wake.dat", orbit_length, bunch_spacing,zgrid, "circular",10); // four columns file
BOOST_CHECK_CLOSE(imped.get_z_coord()[0], 0.00614575, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_coord()[1], 0.0136406, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_coord()[2], 0.0241333, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_coord()[3], 0.037624, tolerance);

BOOST_CHECK_CLOSE(imped.get_z_wake()[0],  9.30907e+08, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_wake()[1],  5.79371e+08, tolerance);
BOOST_CHECK_CLOSE(imped.get_z_wake()[2],  3.78051e+08 , tolerance);
BOOST_CHECK_CLOSE(imped.get_z_wake()[3],  2.55723e+08, tolerance);

BOOST_CHECK_CLOSE(imped.get_x_wake()[0],   -7.40775e+10   , tolerance);
BOOST_CHECK_CLOSE(imped.get_x_wake()[1],   -8.2791e+10, tolerance);
BOOST_CHECK_CLOSE(imped.get_x_wake()[2],   -8.95615e+10 , tolerance);
BOOST_CHECK_CLOSE(imped.get_x_wake()[3],   -9.39278e+10  , tolerance);

BOOST_CHECK_CLOSE(imped.get_y_wake()[0],   -1.26537e+11   , tolerance);
BOOST_CHECK_CLOSE(imped.get_y_wake()[1],   -1.4772e+11 , tolerance);
BOOST_CHECK_CLOSE(imped.get_y_wake()[2],   -1.64615e+11  , tolerance);
BOOST_CHECK_CLOSE(imped.get_y_wake()[3],   -1.77089e+11  , tolerance);

}

BOOST_FIXTURE_TEST_CASE(test_apply, Fixture)
  {
    Bunch bunch(reference_particle, total_num, real_num, comm);
    dummy_populate(bunch);
    Step step(step_length);
    double time_step=10.;
    Impedance imped("test_wake.dat", orbit_length, bunch_spacing, zgrid, "circular",10); // four columns file
    imped.apply(bunch, time_step, step);
}
