#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_track.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-11;

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 9;
const double real_num = 2.0e12;
const int turns = 17;
const double turn_length = 246.8;
const double partial_s = 123.4;

void
dummy_populate_single(Bunch &bunch)
{
    for (int i = 0; i < 6; i += 1) {
        bunch.get_local_particles()[0][i] = 0;
    }
    bunch.get_local_particles()[0][Bunch::id] = 0;
}

struct Fixture
{
    Fixture() :
            reference_particle(pconstants::electron_charge, mass,
                    total_energy), comm_sptr(new Commxx),
        bunch_sptr(new Bunch(reference_particle, 1, real_num, comm_sptr))
    {
        BOOST_TEST_MESSAGE("setup fixture");
        dummy_populate_single(*bunch_sptr);
        bunch_sptr->get_reference_particle().set_trajectory(turns, turn_length,
                0.0);
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch_sptr bunch_sptr;
};

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
    Diagnostics_track diagnostics("dummy.h5", 0);
}

BOOST_FIXTURE_TEST_CASE(is_serial, Fixture)
{
    Diagnostics_track diagnostics("dummy.h5", 0);
    BOOST_CHECK(diagnostics.is_serial());
}


BOOST_FIXTURE_TEST_CASE(write_, Fixture)
{
    Diagnostics_track diagnostics("dummy.h5", 0);
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();
}

BOOST_FIXTURE_TEST_CASE(write_track_sin_x, Fixture)
{
    Diagnostics_track diagnostics("dummy.h5", 0);
    diagnostics.set_bunch_sptr(bunch_sptr);
    double length = 0.1;
    for (int i = 0; i < 200; ++i) {
        bunch_sptr->get_reference_particle().increment_trajectory(length);
        bunch_sptr->get_local_particles()[0][Bunch::x] = sin(
                bunch_sptr->get_reference_particle().get_trajectory_length());
        bunch_sptr->get_local_particles()[0][Bunch::xp] = cos(
                bunch_sptr->get_reference_particle().get_trajectory_length());
        diagnostics.update();
        diagnostics.write();
    }
}

BOOST_FIXTURE_TEST_CASE(serialize_, Fixture)
{
    Diagnostics_track diagnostics("dummy.h5", 0);
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();

    xml_save(diagnostics, "diagnostics_track.xml");

    Diagnostics_track loaded;
    xml_load(loaded, "diagnostics_track.xml");
}
