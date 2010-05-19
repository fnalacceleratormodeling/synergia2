#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/diagnostics.h"
#include "components/foundation/physical_constants.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 9;
const double real_num = 2.0e12;
const int turns = 17;
const double turn_length = 246.8;
const double partial_s = 123.4;

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
        reference_particle(constants::electron_charge, mass, total_energy),
                comm(MPI_COMM_WORLD), bunch(reference_particle, total_num,
                        real_num, comm)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        dummy_populate(bunch);
        bunch.get_reference_particle().set_trajectory(turns, turn_length,
                partial_s);
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
};

BOOST_AUTO_TEST_CASE(construct)
{
    Diagnostics
    diagnostics();
}

BOOST_FIXTURE_TEST_CASE(construct2, Fixture)
{
    Diagnostics diagnostics(bunch);
}

BOOST_FIXTURE_TEST_CASE(get_s, Fixture)
{
    Diagnostics diagnostics(bunch);
    BOOST_CHECK_CLOSE(diagnostics.get_s(), partial_s, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_repetition, Fixture)
{
    Diagnostics diagnostics(bunch);
    BOOST_CHECK_EQUAL(diagnostics.get_repetition(), turns);
}

BOOST_FIXTURE_TEST_CASE(get_trajectory_length, Fixture)
{
    Diagnostics diagnostics(bunch);
    BOOST_CHECK_CLOSE(diagnostics.get_trajectory_length(),
            turns * turn_length + partial_s, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_mean, Fixture)
{
    Diagnostics diagnostics(bunch);
#include "test_diagnostics_get_mean.icc"
}

BOOST_FIXTURE_TEST_CASE(get_std, Fixture)
{
    Diagnostics diagnostics(bunch);
#include "test_diagnostics_get_std.icc"
}

BOOST_FIXTURE_TEST_CASE(init_writers, Fixture)
{
    Diagnostics diagnostics(bunch);
    hid_t hdf5_file = H5Fcreate("test_init_writers.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    H5Fclose(hdf5_file);
}

BOOST_FIXTURE_TEST_CASE(write_hdf5_no_init, Fixture)
{
    Diagnostics diagnostics(bunch);
    hid_t hdf5_file = H5Fcreate("test_write_hdf5_no_init.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    bool caught_error = false;
    try {
        diagnostics.write_hdf5();
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    H5Fclose(hdf5_file);
    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(write_hdf5, Fixture)
{
    Diagnostics diagnostics(bunch);
    hid_t hdf5_file = H5Fcreate("test_write_hdf5.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    diagnostics.write_hdf5();
    H5Fclose(hdf5_file);
}
// test_note: We are not (yet) testing the content of the output file.

// n.b. no test for update because it is called internally for other tests.

BOOST_AUTO_TEST_CASE(construct_full2)
{
    Diagnostics_full2 diagnostics;
}

BOOST_FIXTURE_TEST_CASE(construct2_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
}

BOOST_FIXTURE_TEST_CASE(get_s_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
    BOOST_CHECK_CLOSE(diagnostics.get_s(), partial_s, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_mean_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
#include "test_diagnostics_get_mean.icc"
}

BOOST_FIXTURE_TEST_CASE(get_std_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
#include "test_diagnostics_get_mean.icc"
}

BOOST_FIXTURE_TEST_CASE(get_mom2_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
#include "test_diagnostics_get_mom2.icc"
}

BOOST_FIXTURE_TEST_CASE(get_corr_full2, Fixture)
{
    const double tolerance_corr = 1.0e-10;
    Bunch bunch2(reference_particle, total_num, real_num, comm);
    MArray2d_ref particles(bunch2.get_local_particles());
#include "test_diagnostics_get_random_particles.icc"
    Diagnostics_full2 diagnostics(bunch2);
#include "test_diagnostics_get_corr.icc"
}

const double tolerance_emit2d = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(get_emitx_full2, Fixture)
{
    Bunch bunch2(reference_particle, total_num, real_num, comm);
    MArray2d_ref particles(bunch2.get_local_particles());
#include "test_diagnostics_get_random_particles.icc"
    Diagnostics_full2 diagnostics(bunch2);
#include "test_diagnostics_get_emitx.icc"
}

BOOST_FIXTURE_TEST_CASE(get_emity_full2, Fixture)
{
    Bunch bunch2(reference_particle, total_num, real_num, comm);
    MArray2d_ref particles(bunch2.get_local_particles());
#include "test_diagnostics_get_random_particles.icc"
    Diagnostics_full2 diagnostics(bunch2);
#include "test_diagnostics_get_emity.icc"
}

BOOST_FIXTURE_TEST_CASE(get_emitz_full2, Fixture)
{
    Bunch bunch2(reference_particle, total_num, real_num, comm);
    MArray2d_ref particles(bunch2.get_local_particles());
#include "test_diagnostics_get_random_particles.icc"
    Diagnostics_full2 diagnostics(bunch2);
#include "test_diagnostics_get_emitz.icc"
}

const double tolerance_emit4d = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(get_emitxy_full2, Fixture)
{
    Bunch bunch2(reference_particle, total_num, real_num, comm);
    MArray2d_ref particles(bunch2.get_local_particles());
#include "test_diagnostics_get_random_particles.icc"
    Diagnostics_full2 diagnostics(bunch2);
#include "test_diagnostics_get_emitxy.icc"
}

const double tolerance_emit6d = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(get_emitxyz_full2, Fixture)
{
    Bunch bunch2(reference_particle, total_num, real_num, comm);
    MArray2d_ref particles(bunch2.get_local_particles());
#include "test_diagnostics_get_random_particles.icc"
    Diagnostics_full2 diagnostics(bunch2);
#include "test_diagnostics_get_emitxyz.icc"
}

BOOST_FIXTURE_TEST_CASE(init_writers_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
    hid_t hdf5_file = H5Fcreate("test_init_writers_full2.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    H5Fclose(hdf5_file);
}

BOOST_FIXTURE_TEST_CASE(write_hdf5_full2, Fixture)
{
    Diagnostics_full2 diagnostics(bunch);
    hid_t hdf5_file = H5Fcreate("test_write_hdf5_full2.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    diagnostics.init_writers(hdf5_file);
    diagnostics.write_hdf5();
    H5Fclose(hdf5_file);
}
// test_note: We are not (yet) testing the content of the output file.
