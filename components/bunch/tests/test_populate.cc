#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/bunch/populate.h"
#include "components/bunch/diagnostics.h"
#include "utils/boost_test_mpi_fixture.h"
#include "utils/multi_array_typedefs.h"
#include "utils/multi_array_print.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

const double mass = 100.0;
const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 1000;
const double real_num = 2.0e12;
const double default_s = 123.4;

void
compare_multi_array(Const_MArray1d_ref const& a, Const_MArray1d_ref const& b,
        double tolerance)
{
    BOOST_CHECK_EQUAL(a.shape()[0],b.shape()[0]);
    for (int i = 0; i < a.shape()[0]; ++i) {
        if (a[i] == 0.0) {
            BOOST_CHECK_SMALL(b[i], tolerance);
        } else {
            BOOST_CHECK_CLOSE(a[i], b[i], tolerance);
        }
    }
}

void
compare_multi_array(Const_MArray2d_ref const& a, Const_MArray2d_ref const& b,
        double tolerance)
{
    BOOST_CHECK_EQUAL(a.shape()[0],b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1],b.shape()[1]);
    for (int i = 0; i < a.shape()[0]; ++i) {
        for (int j = 0; j < a.shape()[1]; ++j) {
            if (a[i][j] == 0.0) {
                BOOST_CHECK_SMALL(b[i][j], tolerance);
            } else {
                BOOST_CHECK_CLOSE(a[i][j], b[i][j], tolerance);
            }
        }
    }
}

struct Fixture
{
    Fixture() :
        reference_particle(mass, total_energy), comm(MPI_COMM_WORLD), bunch(
                reference_particle, proton_charge, total_num, real_num, comm),
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
    Diagnostics_full2 diagnostics(bunch, s);
    compare_multi_array(means, diagnostics.get_mean(), tolerance);
    compare_multi_array(covariances, diagnostics.get_mom2(), tolerance);
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
    }
    populate_6d(distribution, bunch, means, covariances);
    Diagnostics_full2 diagnostics(bunch, s);
    compare_multi_array(means, diagnostics.get_mean(), tolerance);
    compare_multi_array(covariances, diagnostics.get_mom2(), tolerance);
}
