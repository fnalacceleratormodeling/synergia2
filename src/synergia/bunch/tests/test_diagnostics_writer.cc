#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics_writer.h"
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
        reference_particle(pconstants::electron_charge, mass, total_energy),
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

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_construct.h5",
            diagnostics_sptr);
}

BOOST_AUTO_TEST_CASE(construct_dummy)
{
    Diagnostics_writer dummy_diagnostics;
}

BOOST_FIXTURE_TEST_CASE(is_dummy_false, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_is_dummy_false.h5",
            diagnostics_sptr);
    BOOST_CHECK(!diagnostics_writer.is_dummy());
}

BOOST_AUTO_TEST_CASE(is_dummy_true)
{
    Diagnostics_writer diagnostics_writer;
    BOOST_CHECK(diagnostics_writer.is_dummy());
}

BOOST_FIXTURE_TEST_CASE(construct_full2, Fixture)
{
    Diagnostics_full2_sptr diagnostics_sptr(new Diagnostics_full2(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_construct_full2.h5",
            diagnostics_sptr);
}

BOOST_FIXTURE_TEST_CASE(get_diagnostics_sptr, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer(
            "test_writer_get_diagnostics_sptr.h5", diagnostics_sptr);
    Diagnostics_sptr retrieved_d_s = diagnostics_writer.get_diagnostics_sptr();
    BOOST_CHECK_EQUAL(diagnostics_sptr, retrieved_d_s);
}

BOOST_FIXTURE_TEST_CASE(get_count, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_get_count.h5",
            diagnostics_sptr);
    BOOST_CHECK_EQUAL(diagnostics_writer.get_count(), 0);
    diagnostics_writer.update_and_write(bunch);
    diagnostics_writer.update_and_write(bunch);
    diagnostics_writer.update_and_write(bunch);
    BOOST_CHECK_EQUAL(diagnostics_writer.get_count(), 3);
}

BOOST_FIXTURE_TEST_CASE(set_count, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_get_count.h5",
            diagnostics_sptr);
    diagnostics_writer.set_count(17);
    BOOST_CHECK_EQUAL(diagnostics_writer.get_count(), 17);
}

BOOST_FIXTURE_TEST_CASE(write_, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_write.h5",
            diagnostics_sptr);
    diagnostics_writer.get_diagnostics_sptr()->update(bunch);
    diagnostics_writer.write();
}

BOOST_FIXTURE_TEST_CASE(update_and_write, Fixture)
{
    Diagnostics_basic_sptr diagnostics_sptr(new Diagnostics_basic(bunch));
    Diagnostics_writer diagnostics_writer("test_writer_update_and_write.h5",
            diagnostics_sptr);
    diagnostics_writer.update_and_write(bunch);
}

BOOST_FIXTURE_TEST_CASE(update_and_write_full2, Fixture)
{
    Diagnostics_full2_sptr diagnostics_sptr(new Diagnostics_full2(bunch));
    Diagnostics_writer diagnostics_writer(
            "test_writer_update_and_write_full2.h5", diagnostics_sptr);
    diagnostics_writer.update_and_write(bunch);
}

BOOST_FIXTURE_TEST_CASE(update_and_write_particles, Fixture)
{
    Diagnostics_particles_sptr diagnostics_sptr(
            new Diagnostics_particles(bunch));
    Diagnostics_writer diagnostics_writer(
            "test_writer_update_and_write_particles.h5", diagnostics_sptr);
    diagnostics_writer.update_and_write(bunch);
    diagnostics_writer.update_and_write(bunch);
    diagnostics_writer.update_and_write(bunch);
}

BOOST_FIXTURE_TEST_CASE(dummy_do_nothing, Fixture)
{
    Diagnostics_writer dummy;
    dummy.get_diagnostics_sptr();
    dummy.write();
    dummy.update_and_write(bunch);
}

BOOST_AUTO_TEST_CASE(no_diagnostics_)
{
    BOOST_CHECK(no_diagnostics().is_dummy());
}
