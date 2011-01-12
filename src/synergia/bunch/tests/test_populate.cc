#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/multi_array_check_equal.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

const double mass = 100.0;
const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 1000;
const double real_num = 2.0e12;
const double default_s = 123.4;


struct Fixture
{
    Fixture() :
        reference_particle(proton_charge, mass, total_energy), comm(MPI_COMM_WORLD), bunch(
                reference_particle, total_num, real_num, comm),
                distribution(0, comm), s(default_s)
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
    Random_distribution distribution;
    double s;
};

BOOST_FIXTURE_TEST_CASE(populate_6d_diagonal, Fixture)
{
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] = i * 0.1;
        covariances[i][i] = (i + 1) * (i + 1);
    }
    populate_6d(distribution, bunch, means, covariances);
    Diagnostics_full2 diagnostics(bunch);
    multi_array_check_equal(means, diagnostics.get_mean(), tolerance);
    multi_array_check_equal(covariances, diagnostics.get_mom2(), tolerance);
}

BOOST_FIXTURE_TEST_CASE(populate_6d_general, Fixture)
{
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] = i * 7.2;
        for (int j = i; j < 6; ++j) {
            covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1);
        }
        covariances[i][i] *= 10.0; // this makes for a positive-definite matrix
    }
    populate_6d(distribution, bunch, means, covariances);
    Diagnostics_full2 diagnostics(bunch);
    multi_array_check_equal(means, diagnostics.get_mean(), tolerance);
    multi_array_check_equal(covariances, diagnostics.get_mom2(), tolerance);
}
