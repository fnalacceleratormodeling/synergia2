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

BOOST_AUTO_TEST_CASE(construct)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
}

BOOST_AUTO_TEST_CASE(get_particle_charge)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    BOOST_CHECK_CLOSE(bunch.get_particle_charge(),proton_charge,tolerance);
}

BOOST_AUTO_TEST_CASE(set_particle_charge)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    const int electron_charge = -1;
    bunch.set_particle_charge(electron_charge);
    BOOST_CHECK_CLOSE(bunch.get_particle_charge(),electron_charge,tolerance);
}

BOOST_AUTO_TEST_CASE(get_real_num)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    BOOST_CHECK_CLOSE(bunch.get_real_num(),real_num,tolerance);
}

BOOST_AUTO_TEST_CASE(set_real_num)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    const double new_real_num = real_num * 1.5;
    bunch.set_real_num(new_real_num);
    BOOST_CHECK_CLOSE(bunch.get_real_num(),new_real_num,tolerance);
}

BOOST_AUTO_TEST_CASE(get_local_num)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    // n.b.: this test assumes that we are running on one processor
    BOOST_CHECK_EQUAL(bunch.get_local_num(),total_num);
}

BOOST_AUTO_TEST_CASE(set_local_num)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    const int new_local_num = 47;
    bunch.set_local_num(new_local_num);
    BOOST_CHECK_EQUAL(bunch.get_local_num(),new_local_num);
}

BOOST_AUTO_TEST_CASE(get_total_num)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    BOOST_CHECK_EQUAL(bunch.get_total_num(),total_num);
}

BOOST_AUTO_TEST_CASE(update_total_num)
{
    Reference_particle reference_particle(total_energy);
    Bunch bunch(reference_particle, proton_charge, total_num, real_num, Commxx(
            MPI_COMM_WORLD));
    const int new_local_num = 47;
    bunch.set_local_num(new_local_num);
    bunch.update_total_num();
    BOOST_CHECK_EQUAL(bunch.get_total_num(),
            new_local_num*Commxx(MPI_COMM_WORLD).get_size());
}

