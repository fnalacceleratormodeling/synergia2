#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/floating_point.h"

const double tolerance = 1.0e-12;
const double small = 1.0e-14;
const double medium = 1.0;
const double large = 1.0e20;

const std::complex<double > tolerance_complex = std::complex<double >(1.0e-12, 
        1.0e-12);
const std::complex<double > small_complex = std::complex<double >(1.0e-14, 
        1.0e-14); 
const std::complex<double > medium_complex = std::complex<double >(1.0, 1.0);
const std::complex<double > large_complex = std::complex<double >(1.0e20, 
        1.0e20);

BOOST_AUTO_TEST_CASE(equal_medium)
{
    BOOST_CHECK(floating_point_equal(medium, medium+0.5*tolerance, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_large)
{
    BOOST_CHECK(floating_point_equal(large, large + 0.5*large*tolerance, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_small)
{
    BOOST_CHECK(floating_point_equal(small, small*0.5, tolerance));
}

BOOST_AUTO_TEST_CASE(not_equal_medium)
{
    BOOST_CHECK(!floating_point_equal(medium, medium+2.0*tolerance, tolerance));
}

BOOST_AUTO_TEST_CASE(not_equal_large)
{
    BOOST_CHECK(!floating_point_equal(large, large + 2.0*large*tolerance, tolerance));
}

BOOST_AUTO_TEST_CASE(not_equal_small)
{
    BOOST_CHECK(!floating_point_equal(small, small+2.0*tolerance, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_medium_negative)
{
    BOOST_CHECK(floating_point_equal(-medium, -medium+0.5*tolerance, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_medium_complex1)
{
    BOOST_CHECK(complex_floating_point_equal(medium_complex, medium_complex+std::complex<double >(0.5, 0.0)*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_medium_complex2)
{
    BOOST_CHECK(complex_floating_point_equal(medium_complex, medium_complex+std::complex<double >(0.0, 0.5)*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_large_complex1)
{
    BOOST_CHECK(complex_floating_point_equal(large_complex, large_complex+std::complex<double >(0.5, 0.0)*large_complex*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_large_complex2)
{
    BOOST_CHECK(complex_floating_point_equal(large_complex, large_complex+std::complex<double >(0.0, 0.5)*large_complex*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_small_complex1)
{
    BOOST_CHECK(complex_floating_point_equal(small_complex, small_complex*std::complex<double >(0.5, 0.0), tolerance));
}

BOOST_AUTO_TEST_CASE(equal_small_complex2)
{
    BOOST_CHECK(complex_floating_point_equal(small_complex, small_complex*std::complex<double >(0.0, 0.5), tolerance));
}

BOOST_AUTO_TEST_CASE(not_equal_medium_complex)
{
    BOOST_CHECK(!complex_floating_point_equal(medium_complex, medium_complex+std::complex<double >(2.0, 2.0)*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(not_equal_large_complex)
{
    BOOST_CHECK(!complex_floating_point_equal(large_complex, large_complex+std::complex<double >(2.0, 2.0)*large_complex*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(not_equal_small_complex)
{
    BOOST_CHECK(!complex_floating_point_equal(small_complex, small_complex+std::complex<double >(2.0, 2.0)*tolerance_complex, tolerance));
}

BOOST_AUTO_TEST_CASE(equal_medium_complex_negative)
{
    BOOST_CHECK(complex_floating_point_equal(-medium_complex, -medium_complex+std::complex<double >(0.5, 0.5)*tolerance_complex, tolerance));
}
