#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/complex_error_function.h"

#include <complex>
#include <cmath>
#include <iostream>

const double tolerance = 1.0e-10;
const std::complex<double > z = std::complex<double >(1.0, 1.0);

BOOST_AUTO_TEST_CASE(check_overflow)
{
    std::complex<double > w;
    std::complex<double > a = std::complex<double > (1.0, 1.0);
    BOOST_CHECK(wofz(a).real());
    BOOST_CHECK(wofz(a).imag());
    double large_number = 1.0e154;
    a = std::complex<double > (large_number, 0.0);
    bool caught_error(false);
    try {
        wofz(a);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);

    a = std::complex<double > (0.0, large_number);
    caught_error = false;
    try {
        wofz(a);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);

    a = std::complex<double > (large_number, large_number);
    caught_error = false;
    try {
        wofz(a);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    BOOST_CHECK(caught_error == true);
}
