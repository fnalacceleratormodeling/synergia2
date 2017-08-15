#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/hdf5_file.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-11;

const double mass = 100.0;
const double total_energy = 125.0;
const int proton_charge = 1;
const int total_num = 10000;
const double real_num = 2.0e12;
const double default_s = 123.4;
const int seed = 987654321;

struct Fixture
{
    Fixture() :
        reference_particle(proton_charge, mass, total_energy), comm_sptr(
                new Commxx), bunch(reference_particle, total_num, real_num,
                comm_sptr), distribution(seed, *comm_sptr), s(default_s)
    {
        BOOST_TEST_MESSAGE("setup fixture");
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Reference_particle reference_particle;
    Commxx_sptr comm_sptr;
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

   covariances[1][1] *=0.00001;
   covariances[3][3] *=0.00001;
   covariances[5][5]*=0.00001;




    populate_6d(distribution, bunch, means, covariances);
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    multi_array_check_equal(means, bunch_mean, tolerance);
    multi_array_check_equal(covariances, bunch_mom2, tolerance);
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
    for (int i = 0; i < 6; ++i) {
         covariances[1][i] *=0.001;
         covariances[i][1] *=0.001;
         covariances[3][i] *=0.001;
         covariances[i][3] *=0.001;
         covariances[5][i] *=0.001;
         covariances[i][5] *=0.001;
    }



    populate_6d(distribution, bunch, means, covariances);
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    multi_array_check_equal(means, bunch_mean, tolerance);
    multi_array_check_equal(covariances, bunch_mom2, tolerance);
}

BOOST_FIXTURE_TEST_CASE(populate_6d_bad_shapes, Fixture)
{
    MArray2d covariances(boost::extents[6][6]);
    MArray2d covariances_bad(boost::extents[6][7]);
    MArray1d means(boost::extents[6]);
    MArray1d means_bad(boost::extents[5]);

    bool caught = false;
    try {
        populate_6d(distribution, bunch, means_bad, covariances);
    }
    catch (std::runtime_error &) {
        caught = true;
    }
    BOOST_CHECK(caught);

    caught = false;
    try {
        populate_6d(distribution, bunch, means, covariances_bad);
    }
    catch (std::runtime_error &) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_FIXTURE_TEST_CASE(populate_6d_truncated_diagonal, Fixture)
{
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    MArray1d limits(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] = i * 0.1;
        covariances[i][i] = (i + 1) * (i + 1);
        limits[i] = 3.0;
    }

    covariances[1][1] *= 0.00001;
    covariances[3][3] *= 0.00001;
    covariances[5][5] *= 0.00001;

    populate_6d_truncated(distribution, bunch, means, covariances, limits);
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    multi_array_check_equal(means, bunch_mean, tolerance);
    multi_array_check_equal(covariances, bunch_mom2, tolerance);

    bool found_miscreant(false);
    for (int particle = 0; particle < bunch.get_local_num(); ++particle) {
        for (int i = 0; i < 6; ++i) {
            if (std::abs(bunch.get_local_particles()[particle][i]-means[i])
                    > std::sqrt(covariances[i][i]) * limits[i]) {
                found_miscreant = true;
            }
        }
    }
    BOOST_CHECK(!found_miscreant);
}

BOOST_FIXTURE_TEST_CASE(populate_6d_truncated_general, Fixture)
{
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    MArray1d limits(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] = i * 7.2;
        for (int j = i; j < 6; ++j) {
            covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1);
        }
        covariances[i][i] *= 10.0; // this makes for a positive-definite matrix
        limits[i] = 3.0;
    }
    for (int i = 0; i < 6; ++i) {
         covariances[1][i] *=0.001;
         covariances[i][1] *=0.001;
         covariances[3][i] *=0.001;
         covariances[i][3] *=0.001;
         covariances[5][i] *=0.001;
         covariances[i][5] *=0.001;
    }

    populate_6d_truncated(distribution, bunch, means, covariances, limits);
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    multi_array_check_equal(means, bunch_mean, tolerance);
    multi_array_check_equal(covariances, bunch_mom2, tolerance);
}

BOOST_FIXTURE_TEST_CASE(populate_transverse_gaussian_general, Fixture)
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
    const double cdt = 3.14;
    populate_transverse_gaussian(distribution, bunch, means, covariances, cdt);
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Core_diagnostics::calculate_mom2(bunch, bunch_mean));
    for (int i = 0; i < 4; ++i) {
        if (std::abs(means[i]) < 10 * tolerance) {
            BOOST_CHECK_SMALL(bunch_mean[i], tolerance);
        } else {
            BOOST_CHECK_CLOSE(means[i], bunch_mean[i], tolerance);
        }
        for (int j = i; j < 4; ++j) {
            BOOST_CHECK_CLOSE(covariances[i][j], bunch_mom2[i][j],
                    tolerance);
        }
    }
    for (int i = 0; i < 6; ++i) {
        if (i != 4) {
            BOOST_CHECK_SMALL(bunch_mom2[i][4], tolerance);
            BOOST_CHECK_CLOSE(covariances[i][5], bunch_mom2[i][5],
                    tolerance);
            BOOST_CHECK_CLOSE(covariances[5][i], bunch_mom2[5][i],
                    tolerance);
        }
    }
    BOOST_CHECK_SMALL(bunch_mean[4], tolerance);
    BOOST_CHECK_CLOSE(means[5], bunch_mean[5], tolerance);
    BOOST_CHECK_CLOSE(cdt*cdt/12.0, bunch_mom2[4][4],
            tolerance);

    double sum2 = 0.0;
    double sum3 = 0.0;
    double sum4 = 0.0;
    double sum5 = 0.0;
    double sum6 = 0.0;
    MArray2d_ref particles(bunch.get_local_particles());
    Hdf5_file f("particles.h5", Hdf5_file::truncate);
    f.write(particles, "particles");
    int num = bunch.get_local_num();
    for (int p = 0; p < num; ++p) {
        double cdt = particles[p][Bunch::cdt]
                - bunch_mean[Bunch::cdt];
        double cdt2 = cdt * cdt;
        sum2 += cdt2;
        double cdt3 = cdt2 * cdt;
        sum3 += cdt3;
        double cdt4 = cdt3 * cdt;
        sum4 += cdt4;
        double cdt5 = cdt4 * cdt;
        sum5 += cdt5;
        double cdt6 = cdt5 * cdt;
        sum6 += cdt6;
    }
    BOOST_CHECK_CLOSE(sum2/num, cdt*cdt/12.0, tolerance);
    // n.b. the tolerances for the higher moments are very large
    BOOST_CHECK_SMALL(sum3/num, 0.03);
    BOOST_CHECK_CLOSE(sum4/num, cdt*cdt*cdt*cdt/80.0, 1.0);
    BOOST_CHECK_SMALL(sum5/num, 0.1);
    BOOST_CHECK_CLOSE(sum6/num, cdt*cdt*cdt*cdt*cdt*cdt/448.0, 5.0);
}

BOOST_FIXTURE_TEST_CASE(populate_uniform_cylinder_general, Fixture)
{
    const double radius = 2.718;
    const double cdt = 3.14;
    const double stdxp = 1.2e-3;
    const double stdyp = 3.4e-3;
    const double stddpop = 5.6e-4;
    populate_uniform_cylinder(distribution, bunch, radius, cdt, stdxp, stdyp,
            stddpop);
    MArray1d bunch_mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d bunch_std(Core_diagnostics::calculate_std(bunch, bunch_mean));
    // yes, this tolerance is very large. There are no corrections for
    // finite statistics in populate_uniform_cylinder.
    const double tolerance = 2.0;
    for (int i = 0; i < 6; ++i) {
        BOOST_CHECK_SMALL(bunch_mean[i], tolerance);
    }
    BOOST_CHECK_CLOSE(bunch_std[Bunch::xp], stdxp, tolerance);
    BOOST_CHECK_CLOSE(bunch_std[Bunch::yp], stdyp, tolerance);
    BOOST_CHECK_CLOSE(bunch_std[Bunch::dpop], stddpop, tolerance);

    double sum = 0.0;
    double max = 0.0;
    MArray2d_ref particles(bunch.get_local_particles());
    Hdf5_file f("cylparticles.h5", Hdf5_file::truncate);
    f.write(particles, "particles");
    int num = bunch.get_local_num();
    for (int p = 0; p < num; ++p) {
        double r = std::sqrt(particles[p][Bunch::x] * particles[p][Bunch::x]
                + particles[p][Bunch::y] * particles[p][Bunch::y]);
        if (r>max) {
            max = r;
        }
        sum += r;
    }
    BOOST_CHECK_CLOSE(max, radius, tolerance);
    BOOST_CHECK_CLOSE(sum/num, 2.0*radius/3.0, tolerance);
}

///==================================================================

// populate a KV distribution with specific emittances and lattice parameters.
// verify that the resulting distribution reflects those values.
BOOST_FIXTURE_TEST_CASE(populate_KVlongitudinal, Fixture)
{
    const double kvtolerance = 1.5;  // tolerance for populate tests is large because
    // of finite statistics.  Also notice that alphax and alphay are not zero because
    // statistics of the distributions will not get you correlations close to zero.
    const double emitx = 5.0e-6; // 5 mm-mrad
    const double betax = 10.0;
    const double alphax = 0.75;
    const double emity = 7.0e-6; // 7 mm-mrad
    const double betay = 60.0;
    const double alphay = -1.2;

    populate_transverse_KV_GaussLong(distribution, bunch, emitx, alphax, betax,
                                     emity, alphay, betay, 0.1, 0.0);

    MArray1d mean(Core_diagnostics::calculate_mean(bunch));
    MArray1d std(Core_diagnostics::calculate_std(bunch, mean));
    MArray2d mom2(Core_diagnostics::calculate_mom2(bunch, mean));

    double temitx = std::sqrt(mom2[0][0]*mom2[1][1]-mom2[0][1]*mom2[1][0]);
    double temity = std::sqrt(mom2[2][2]*mom2[3][3]-mom2[2][3]*mom2[3][2]);
    double tmomxx = mom2[0][0];
    double tmomxxp = mom2[0][1];
    double tmomxpxp = mom2[1][1];
    double tmomyy = mom2[2][2];
    double tmomyyp = mom2[2][3];
    double tmomypyp = mom2[3][3];

    BOOST_CHECK_CLOSE(tmomxx, emitx*betax, kvtolerance);
    BOOST_CHECK_CLOSE(tmomxxp, -alphax*emitx, kvtolerance);
    BOOST_CHECK_CLOSE(tmomxpxp, (1.0+alphax*alphax)*emitx/betax, kvtolerance);
    BOOST_CHECK_CLOSE(temitx, emitx, kvtolerance);

    BOOST_CHECK_CLOSE(tmomyy, emity*betay, kvtolerance);
    BOOST_CHECK_CLOSE(tmomyyp, -alphay*emity, kvtolerance);
    BOOST_CHECK_CLOSE(tmomypyp, (1.0+alphay*alphay)*emity/betay, kvtolerance);
    BOOST_CHECK_CLOSE(temity, emity, kvtolerance);

}
