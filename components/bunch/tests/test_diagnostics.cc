// python generated output for this file comes from
// test_diagnostics_crosscheck.py
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/diagnostics.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

const double mass = 100.0;
const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 9;
const double real_num = 2.0e12;
const double default_s = 123.4;

void
dummy_populate(Bunch &bunch)
{
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; i += 1) {
            bunch.get_local_particles()[part][i] = 10.0 * part + 1.1 * i;
        }
        bunch.get_local_particles()[part][Bunch::id] = part;
    }
}

struct Fixture
{
    Fixture() :
        reference_particle(mass, total_energy), comm(MPI_COMM_WORLD), bunch(
                reference_particle, proton_charge, total_num, real_num, comm),
                s(default_s)
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
};

BOOST_AUTO_TEST_CASE(construct)
{
    Diagnostics
    diagnostics();
}

BOOST_FIXTURE_TEST_CASE(construct2, Fixture)
{
    Diagnostics diagnostics(bunch, s);
}

BOOST_FIXTURE_TEST_CASE(get_s, Fixture)
{
    Diagnostics diagnostics(bunch, s);
    BOOST_CHECK_CLOSE(diagnostics.get_s(),default_s,tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_mean, Fixture)
{
    Diagnostics diagnostics(bunch, s);
    // python output for get_mean
    BOOST_CHECK_CLOSE(diagnostics.get_mean()[0], 40.00000, tolerance);
    BOOST_CHECK_CLOSE(diagnostics.get_mean()[1], 41.10000, tolerance);
    BOOST_CHECK_CLOSE(diagnostics.get_mean()[2], 42.20000, tolerance);
    BOOST_CHECK_CLOSE(diagnostics.get_mean()[3], 43.30000, tolerance);
    BOOST_CHECK_CLOSE(diagnostics.get_mean()[4], 44.40000, tolerance);
    BOOST_CHECK_CLOSE(diagnostics.get_mean()[5], 45.50000, tolerance);
    // end python output for get_mean
}

