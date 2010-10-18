#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/floating_point.h"

const double tolerance = 1.0e-12;
const double small = 1.0e-14;
const double medium = 1.0;
const double large = 1.0e20;

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

