#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/diagnostics.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-15;

const double mass = 100.0;
const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 100;
const double real_num = 2.0e12;
const double default_s = 123.4;

void
dummy_populate(Bunch &bunch, int id_offset = 0)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        int id = part + id_offset;
        // coordinates
        for (int i = 0; i < 6; i += 2) {
            bunch.get_local_particles()[part][i] = 10.0 * id + i;
        }
        // momenta
        for (int i = 1; i < 6; i += 2) {
            bunch.get_local_particles()[part][i] = 1e-4 * (10.0 * id + i);
        }
        bunch.get_local_particles()[part][Bunch::id] = id;
    }
}
struct Fixture
{
    Fixture() :
        reference_particle(mass, total_energy), comm(MPI_COMM_WORLD), bunch(
                reference_particle, proton_charge, total_num, real_num, comm),
                s(default_s), diagnostics(bunch,s)
    {
        BOOST_TEST_MESSAGE("setup fixture");
        dummy_populate(bunch);
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Reference_particle reference_particle;
    Commxx comm;
    Bunch bunch;
    double s;
    Diagnostics diagnostics;
};

BOOST_AUTO_TEST_CASE(construct)
{
    Diagnostics diagnostics();
}

BOOST_FIXTURE_TEST_CASE(construct2, Fixture)
{
}
