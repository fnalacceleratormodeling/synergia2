#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/bunch.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-15;

const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 100;
const double real_num = 2.0e12;

struct Fixture
{
    Fixture() :
        reference_particle(total_energy), comm(MPI_COMM_WORLD), bunch(
                reference_particle, proton_charge, total_num, real_num, comm)
    {
        BOOST_TEST_MESSAGE("setup fixture");
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
}

BOOST_FIXTURE_TEST_CASE(get_particle_charge, Fixture)
{
    BOOST_CHECK_CLOSE(bunch.get_particle_charge(),proton_charge,tolerance);
}

BOOST_FIXTURE_TEST_CASE(set_particle_charge, Fixture)
{
    const int electron_charge = -1;
    bunch.set_particle_charge(electron_charge);
    BOOST_CHECK_CLOSE(bunch.get_particle_charge(),electron_charge,tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_real_num, Fixture)
{
    BOOST_CHECK_CLOSE(bunch.get_real_num(),real_num,tolerance);
}

BOOST_FIXTURE_TEST_CASE(set_real_num, Fixture)
{
    const double new_real_num = real_num * 1.5;
    bunch.set_real_num(new_real_num);
    BOOST_CHECK_CLOSE(bunch.get_real_num(),new_real_num,tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_local_num, Fixture)
{
    // n.b.: this test assumes that we are running on one processor
    BOOST_CHECK_EQUAL(bunch.get_local_num(),total_num);
}

BOOST_FIXTURE_TEST_CASE(set_local_num, Fixture)
{
    int new_local_num = total_num / comm.get_size() - 5;
    bunch.set_local_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_num(),new_local_num);
}

BOOST_FIXTURE_TEST_CASE(get_total_num, Fixture)
{
    BOOST_CHECK_EQUAL(bunch.get_total_num(),total_num);
}

BOOST_FIXTURE_TEST_CASE(update_total_num, Fixture)
{
    const int new_local_num = 47;
    bunch.set_local_num(new_local_num);
    bunch.update_total_num();
    BOOST_CHECK_EQUAL(bunch.get_total_num(),
            new_local_num*Commxx(MPI_COMM_WORLD).get_size());
}

BOOST_FIXTURE_TEST_CASE(get_local_particles, Fixture)
{
    MArray2d_ref local_particles(bunch.get_local_particles());
    BOOST_CHECK_EQUAL(local_particles.shape()[1],7);
    BOOST_CHECK(local_particles.shape()[0] >= bunch.get_local_num());
}

BOOST_FIXTURE_TEST_CASE(increase_local_num, Fixture)
{
    const int small_total_num = 10;
    const int increase = 5;
    Bunch bunch2(reference_particle, proton_charge, small_total_num, real_num,
            comm);
    // populate bunch2
    MArray2d_ref particles(bunch2.get_local_particles());
    int old_local_num = bunch2.get_local_num();
    for (int particle = 0; particle < old_local_num; ++particle) {
        for (int index = 0; index < 7; ++index) {
            particles[particle][index] = particle * 10.0 + index;
        }
    }

    // expand bunch2 and verify that old values are still there
    bunch2.set_local_num(old_local_num+increase);
    MArray2d_ref particles2(bunch2.get_local_particles());
    BOOST_CHECK_EQUAL(particles2.shape()[1],7);
    BOOST_CHECK(particles2.shape()[0] >= bunch2.get_local_num());
    for (int particle = 0; particle < old_local_num; ++particle) {
        for (int index = 0; index < 7; ++index) {
            BOOST_CHECK_CLOSE(particles[particle][index],
                    particles[particle][index],tolerance);
        }
    }
}

BOOST_FIXTURE_TEST_CASE(get_state, Fixture)
{
    Bunch::State state;
    state = bunch.get_state();
    BOOST_CHECK_EQUAL(state,Bunch::fixed_z);
}
