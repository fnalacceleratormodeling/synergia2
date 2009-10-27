#define BOOST_TEST_MAIN
#include <vector>
#include <boost/test/unit_test.hpp>
#include "components/foundation/distribution.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-15;
const unsigned long int test_seed = 12345678;
const int array_length = 10000;
const int default_min = -7.1;
const int default_max = 3.25;
const int n_sigma = 4.0; // maximum number of standard deviations for statistical tests
const double pi = 3.1415926535897932384626433832;

void
verify_uniform_distribution(MArray1d_ref const& x, double min, double max,
        int num_bins = 4)
{
    std::vector<double > bins(num_bins);
    for (int i = 0; i < num_bins; ++i) {
        bins.at(i) = 0.0;
    }
    int N = x.size();
    double found_max = min;
    double found_min = max;
    double out_of_bounds = false;
    for (int i = 0; i < N; ++i) {
        if (x[i] < found_min) {
            found_min = x[i];
        }
        if (x[i] > found_max) {
            found_max = x[i];
        }

        int which_bin = static_cast<int > ((x[i] - min) / (max - min)
                * num_bins);
        if (which_bin < 0 || which_bin >= num_bins) {
            out_of_bounds = true;
        } else {
            bins.at(which_bin) += 1.0;
        }
    }
    BOOST_CHECK(!out_of_bounds);
    // Testing min and max may be redundant with the out_of_bounds test,
    // but it generates a clearer error message in case of failure.
    BOOST_CHECK(found_min >= min);
    BOOST_CHECK(found_max < max);
    double chisq = 0.0;
    double expected = (1.0 * N) / num_bins;
    for (int i = 0; i < num_bins; ++i) {
        chisq += (bins.at(i) - expected) * (bins.at(i) - expected) / bins.at(i);
    }
    double chisq_dof = chisq / (num_bins - 1.0);
    BOOST_CHECK(chisq_dof < n_sigma);
}

void
verify_unit_gaussian_distribution(MArray1d_ref const& x)
{
    double sum = 0;
    double sum2 = 0;
    int N = x.size();
    for (int i = 0; i < N; ++i) {
        sum += x[i];
    }
    double mean = sum / N;
    for (int i = 0; i < N; ++i) {
        sum2 += (x[i] - mean) * (x[i] - mean);
    }
    double stddev = std::sqrt(sum2 / N);
    BOOST_CHECK(std::abs(mean) < n_sigma * std::sqrt(1.0/N));
    BOOST_CHECK(std::abs(stddev-1.0) < n_sigma * std::sqrt(2.0 / N));
}

void
verify_unit_disk_distribution(MArray1d_ref const& x, MArray1d_ref const& y)
{
    int N = x.size();
    MArray1d r(boost::extents[N]);
    MArray1d r2(boost::extents[N]);
    MArray1d phi(boost::extents[N]);
    for (int i = 0; i < N; ++i) {
        r2[i] = x[i] * x[i] + y[i] * y[i];
        r[i] = std::sqrt(r2[i]);
        if (r[i] == 0.0) {
            phi[i] = 0.0;
        } else {
            if (x[i] >= 0.0) {
                if (y[i] >= 0.0) {
                    phi[i] = std::asin(y[i] / r[i]);
                } else {
                    phi[i] = 2* pi + std::asin(y[i] / r[i]);
                }
            } else {
                phi[i] = pi - std::asin(y[i] / r[i]);
            }
        }
    }
    verify_uniform_distribution(r2, 0, 1.0);
    verify_uniform_distribution(phi, 0, 2* pi );
}

BOOST_AUTO_TEST_CASE(construct_random)
{
    Random_distribution distribution(0, Commxx(MPI_COMM_WORLD));
}

BOOST_AUTO_TEST_CASE(construct2_random)
{
    Random_distribution distribution(test_seed, Commxx(MPI_COMM_WORLD),
            Random_distribution::ranlxd2);
}

BOOST_AUTO_TEST_CASE(construct3_random)
{
    Random_distribution distribution(test_seed, Commxx(MPI_COMM_WORLD),
            Random_distribution::mt19937);
}

BOOST_AUTO_TEST_CASE(get_seed_random)
{
    Random_distribution distribution(test_seed, Commxx(MPI_COMM_WORLD));
    BOOST_CHECK_EQUAL(test_seed, distribution.get_original_seed());
}

BOOST_AUTO_TEST_CASE(get_default_seed_random)
{
    BOOST_CHECK(0 != Random_distribution::get_default_seed());
}

BOOST_AUTO_TEST_CASE(get_default_seed2_random)
{
    BOOST_CHECK(0 != Random_distribution::get_default_seed("never_find_this_file"));
}

BOOST_AUTO_TEST_CASE(get_random)
{
    // Admittedly, this is a stupid test
    Random_distribution distribution(0, Commxx(MPI_COMM_WORLD));
    BOOST_CHECK(distribution.get() != distribution.get());
}

BOOST_AUTO_TEST_CASE(fill_uniform_random)
{
    Random_distribution distribution(0, Commxx(MPI_COMM_WORLD));
    MArray1d array(boost::extents[array_length]);
    distribution.fill_uniform(array, default_min, default_max);
    verify_uniform_distribution(array, default_min, default_max);
}

BOOST_AUTO_TEST_CASE(fill_unit_gaussian_random)
{
    Random_distribution distribution(0, Commxx(MPI_COMM_WORLD));
    MArray1d array(boost::extents[array_length]);
    distribution.fill_unit_gaussian(array);
    verify_unit_gaussian_distribution(array);
}

BOOST_AUTO_TEST_CASE(fill_unit_disk_random)
{
    Random_distribution distribution(0, Commxx(MPI_COMM_WORLD));
    MArray1d x_array(boost::extents[array_length]);
    MArray1d y_array(boost::extents[array_length]);
    distribution.fill_unit_disk(x_array, y_array);
    verify_unit_disk_distribution(x_array, y_array);
}
