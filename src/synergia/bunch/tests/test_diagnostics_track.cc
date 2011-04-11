#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics.h"
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
        bunch_sptr(new Bunch(reference_particle, 1, real_num, comm)),
                reference_particle(pconstants::electron_charge, mass,
                        total_energy), comm(MPI_COMM_WORLD)
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
    Commxx comm;
    Bunch_sptr bunch_sptr;
};

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
    Diagnostics_track diagnostics(bunch_sptr, "dummy.h5", 0);
}

BOOST_FIXTURE_TEST_CASE(is_serial, Fixture)
{
    Diagnostics_track diagnostics(bunch_sptr, "dummy.h5", 0);
    BOOST_CHECK(diagnostics.is_serial());
}

BOOST_FIXTURE_TEST_CASE(init_writers, Fixture)
{
    Diagnostics_track diagnostics(bunch_sptr, "dummy.h5", 0);
    hid_t hdf5_file = H5Fcreate("test_track_a_00000.h5", H5F_ACC_TRUNC,
            H5P_DEFAULT, H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    H5Fclose(hdf5_file);
}

BOOST_FIXTURE_TEST_CASE(write_hdf5_no_init, Fixture)
{
    Diagnostics_track diagnostics(bunch_sptr, "dummy.h5", 0);
    bool caught_error = false;
    try {
        diagnostics.write_hdf5();
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(write_hdf5, Fixture)
{
    Diagnostics_track diagnostics(bunch_sptr, "dummy.h5", 0);
    hid_t hdf5_file = H5Fcreate("test_track_b_00000.h5", H5F_ACC_TRUNC,
            H5P_DEFAULT, H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    diagnostics.write_hdf5();
    H5Fclose(hdf5_file);
}

BOOST_FIXTURE_TEST_CASE(write_track_sin_x, Fixture)
{
    Diagnostics_track diagnostics(bunch_sptr, "dummy.h5", 0);
    hid_t hdf5_file = H5Fcreate("test_track_c_00000.h5", H5F_ACC_TRUNC,
            H5P_DEFAULT, H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    double length = 0.1;
    for (int i = 0; i < 200; ++i) {
        bunch_sptr->get_reference_particle().increment_trajectory(length);
        bunch_sptr->get_local_particles()[0][Bunch::x] = sin(
                bunch_sptr->get_reference_particle().get_trajectory_length());
        bunch_sptr->get_local_particles()[0][Bunch::xp] = cos(
                bunch_sptr->get_reference_particle().get_trajectory_length());
        diagnostics.update();
        diagnostics.write_hdf5();
    }
    H5Fclose(hdf5_file);
}
