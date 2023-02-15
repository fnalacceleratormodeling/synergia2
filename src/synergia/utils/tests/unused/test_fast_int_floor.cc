#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/fast_int_floor.h"

#include <cmath>

BOOST_AUTO_TEST_CASE(exact_doubles)
{
    for (double x = -10.0; x <= 10.0; x += 0.5) {
        BOOST_CHECK_EQUAL(fast_int_floor(x), static_cast<int > (floor(x)));
    }
}

BOOST_AUTO_TEST_CASE(inexact_doubles)
{
    for (double x = -5.0; x <= 5.0; x += 0.1) {
        BOOST_CHECK_EQUAL(fast_int_floor(x), static_cast<int > (floor(x)));
    }
}
