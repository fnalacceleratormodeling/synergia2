#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 9;
const double real_num = 2.0e12;
const int turns = 17;
const double turn_length = 246.8;
const double partial_s = 123.4;
const int max_particles = 4;

void
dummy_populate(Bunch &bunch)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; i += 1) {
            bunch.get_local_particles()[part][i] = 10.0 * part + (1.0 + part
                    * part / 1000.0) * i;
        }
        bunch.get_local_particles()[part][Bunch::id] = part;
    }
}

struct Fixture
{
    Fixture() :
                bunch_sptr(
                        new Bunch(reference_particle, total_num, real_num, comm)),
                reference_particle(pconstants::electron_charge, mass,
                        total_energy), comm(MPI_COMM_WORLD)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        dummy_populate(*bunch_sptr);
        bunch_sptr->get_reference_particle().set_trajectory(turns, turn_length,
                partial_s);
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Reference_particle reference_particle;
    Commxx comm;
    Bunch_sptr bunch_sptr;
};

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
    Diagnostics_particles diagnostics(bunch_sptr, "dummy.h5");
}

BOOST_FIXTURE_TEST_CASE(construct2, Fixture)
{
    Diagnostics_particles diagnostics(bunch_sptr, "dummy.h5", max_particles);
}

BOOST_FIXTURE_TEST_CASE(is_serial, Fixture)
{
    Diagnostics_particles diagnostics(bunch_sptr, "dummy.h5");
    BOOST_CHECK(!diagnostics.is_serial());
}

BOOST_FIXTURE_TEST_CASE(update, Fixture)
{
    Diagnostics_particles diagnostics(bunch_sptr, "dummy.h5");
    diagnostics.update();
}

BOOST_FIXTURE_TEST_CASE(write_, Fixture)
{
    Diagnostics_particles diagnostics(bunch_sptr,
            "diagnostics_particles_nompi.h5");
    diagnostics.update();
    diagnostics.write();
}

BOOST_FIXTURE_TEST_CASE(write_min_max, Fixture)
{
    Diagnostics_particles diagnostics(bunch_sptr,
            "diagnostics_particles_nompi_35.h5", 3, 5);
    diagnostics.update();
    diagnostics.write();
}
