#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include "synergia/utils/floating_point.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"



BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double mass = 100.0;
const double total_energy = 125.0;
const double real_num = 2.0e12;

struct Fixture
{
    Fixture(int total_num, int total_spectator_num) :
        four_momentum(mass, total_energy), 
        reference_particle(pconstants::proton_charge, four_momentum),
        comm_sptr(new Commxx()), 
        bunch(reference_particle, total_num, total_spectator_num, real_num, comm_sptr),
        total_num(total_num),
        total_spectator_num(total_spectator_num)
    {
        BOOST_TEST_MESSAGE("setup fixture");
    }

    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Four_momentum four_momentum;
    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
    Bunch bunch;

    int total_num;
    int total_spectator_num;
};

struct Fix_p100_s0 : public Fixture
{ Fix_p100_s0() : Fixture(100, 0) { } };

struct Fix_p100_s10 : public Fixture
{ Fix_p100_s10() : Fixture(100, 10) { } };

void
dummy_populate(Bunch &bunch, int offset = 0)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) 
    {
        // coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            bunch.get_local_particles()[part][i] = 10.0 * (part + offset) + i;
        }

        // momenta
        for (int i = 1; i < 6; i += 2) 
        {
            bunch.get_local_particles()[part][i] = 1e-4 * (10.0 * (part + offset) + i);
        }
    }

    // spectator particles
    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        // coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            bunch.get_local_spectator_particles()[part][i] = 10.0 * (part + offset) + i;
        }

        // momenta
        for (int i = 1; i < 6; i += 2) 
        {
            bunch.get_local_spectator_particles()[part][i] = 1e-4 * (10.0 * (part + offset) + i);
        }
    }

}

void
random_populate(Bunch &bunch)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) 
    {
        for (int i = 0; i < 6; i += 1) 
        {
            bunch.get_local_particles()[part][i] = 
                double(std::rand()) / RAND_MAX;
        }

        bunch.get_local_particles()[part][6] = part;
    }

    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        for (int i = 0; i < 6; i += 1) 
        {
            bunch.get_local_spectator_particles()[part][i] = 
                double(std::rand()) / RAND_MAX;
        }

        bunch.get_local_spectator_particles()[part][6] = part;
    }
}

void
compare_bunches(Bunch &bunch1, Bunch &bunch2, double tolerance,
        bool check_state = true, bool check_ids = true, bool check_bucket_index = true)
{
    BOOST_CHECK_EQUAL(
            bunch1.get_reference_particle().get_total_energy(),
            bunch2.get_reference_particle().get_total_energy() );

    BOOST_CHECK_EQUAL(
            bunch1.get_particle_charge(),
            bunch2.get_particle_charge() );

    BOOST_CHECK_CLOSE(bunch1.get_mass(), bunch2.get_mass(), tolerance);
    BOOST_CHECK_CLOSE(bunch1.get_real_num(), bunch2.get_real_num(), tolerance);
    BOOST_CHECK_EQUAL(bunch1.get_local_num(), bunch2.get_local_num());
    BOOST_CHECK_EQUAL(bunch1.get_total_num(), bunch2.get_total_num());

    if (check_bucket_index) 
    {
        BOOST_CHECK_EQUAL(bunch1.is_bucket_index_assigned(), bunch2.is_bucket_index_assigned());
        if (bunch1.is_bucket_index_assigned()){
            BOOST_CHECK_EQUAL(bunch1.get_bucket_index(), bunch2.get_bucket_index());
        }
    }

    BOOST_CHECK_EQUAL(bunch1.get_z_period_length(), bunch2.get_z_period_length());
    BOOST_CHECK_EQUAL(bunch1.get_longitudinal_aperture_length(), bunch2.get_longitudinal_aperture_length());
    BOOST_CHECK_EQUAL(bunch1.is_z_periodic(), bunch2.is_z_periodic());
    BOOST_CHECK_EQUAL(bunch1.has_longitudinal_aperture(),  bunch2.has_longitudinal_aperture());

    if (check_state) 
    {
        BOOST_CHECK_EQUAL(bunch1.get_state(), bunch2.get_state());
    }

    for (int part = 0; part < bunch1.get_local_num(); ++part) 
    {
        // this loop is unrolled in order to give more meaningful error messages
        // i.e., error messages including which component was tested
        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][0],
                bunch2.get_local_particles()[part][0], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][1],
                bunch2.get_local_particles()[part][1], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][2],
                bunch2.get_local_particles()[part][2], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][3],
                bunch2.get_local_particles()[part][3], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][4],
                bunch2.get_local_particles()[part][4], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_particles()[part][5],
                bunch2.get_local_particles()[part][5], tolerance );

        if (check_ids) 
        {
            BOOST_CHECK_CLOSE(
                    bunch1.get_local_particles()[part][6],
                    bunch2.get_local_particles()[part][6], tolerance );
        }
    }

    for (int part = 0; part < bunch1.get_local_spectator_num(); ++part) 
    {
        // this loop is unrolled in order to give more meaningful error messages
        // i.e., error messages including which component was tested
        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][0],
                bunch2.get_local_spectator_particles()[part][0], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][1],
                bunch2.get_local_spectator_particles()[part][1], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][2],
                bunch2.get_local_spectator_particles()[part][2], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][3],
                bunch2.get_local_spectator_particles()[part][3], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][4],
                bunch2.get_local_spectator_particles()[part][4], tolerance );

        BOOST_CHECK_CLOSE(
                bunch1.get_local_spectator_particles()[part][5],
                bunch2.get_local_spectator_particles()[part][5], tolerance );

        if (check_ids) 
        {
            BOOST_CHECK_CLOSE(
                    bunch1.get_local_spectator_particles()[part][6],
                    bunch2.get_local_spectator_particles()[part][6], tolerance );
        }
    }
}


BOOST_FIXTURE_TEST_CASE(construct_no_spectator, Fix_p100_s0)
{
}

BOOST_FIXTURE_TEST_CASE(construct_spectator, Fix_p100_s10)
{
}

BOOST_FIXTURE_TEST_CASE(spectator_num, Fix_p100_s10)
{
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), 10);
    BOOST_CHECK_EQUAL(bunch.get_total_spectator_num(), 10);
}

BOOST_FIXTURE_TEST_CASE(get_local_spectator_particles, Fix_p100_s10)
{
    MArray2d_ref particles(bunch.get_local_spectator_particles());
    BOOST_CHECK_EQUAL(particles.shape()[1], 7);

    unsigned int u_local_num = bunch.get_local_spectator_num();
    BOOST_CHECK(particles.shape()[0] >= u_local_num);

    unsigned int u_local_num_slots = bunch.get_local_spectator_num_slots();
    BOOST_CHECK(particles.shape()[0] == u_local_num_slots);
}

BOOST_FIXTURE_TEST_CASE(get_const_local_spectator_particles, Fix_p100_s10)
{
    Const_MArray2d_ref particles(bunch.get_local_spectator_particles());
    BOOST_CHECK_EQUAL(particles.shape()[1], 7);

    unsigned int u_local_num = bunch.get_local_spectator_num();
    BOOST_CHECK(particles.shape()[0] >= u_local_num);

    unsigned int u_local_num_slots = bunch.get_local_spectator_num_slots();
    BOOST_CHECK(particles.shape()[0] == u_local_num_slots);
}

BOOST_FIXTURE_TEST_CASE(set_local_spectator_num, Fix_p100_s10)
{
    int new_local_num = total_spectator_num / comm_sptr->get_size() - 5;

    bunch.set_local_spectator_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), 100);
}

BOOST_FIXTURE_TEST_CASE(set_local_spectator_num_neg, Fix_p100_s10)
{
    int new_local_num = total_spectator_num / comm_sptr->get_size() - 12;

    bunch.set_local_spectator_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_spectator_num(), 0);
    BOOST_CHECK_EQUAL(bunch.get_local_num(), 100);
}

BOOST_FIXTURE_TEST_CASE(get_total_spectator_num, Fix_p100_s10)
{
    BOOST_CHECK_EQUAL(bunch.get_total_spectator_num(), 10);
}

BOOST_FIXTURE_TEST_CASE(update_total_num_with_new_spectator, Fix_p100_s10)
{
    const int new_local_num = bunch.get_local_num() - 7;
    bunch.set_local_spectator_num(new_local_num);
    bunch.update_total_num();

    BOOST_CHECK_EQUAL(
            bunch.get_total_spectator_num(), 
            new_local_num*Commxx().get_size() );

    BOOST_CHECK_EQUAL(bunch.get_total_num(), 100);
}
