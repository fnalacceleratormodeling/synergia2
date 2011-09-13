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

// csp:Before running this test, enable the part of Leo's code in complex_error_function.h.
BOOST_AUTO_TEST_CASE(check_w_function)
{
    std::complex<double > w1, w2;
    for(int i = 0; i <= 100; ++i) {
        double x = -1.0 + i * 2.0 / 100.0;
        for(int j = 0; j <= 100; ++j) {
            double y = -1.0 + j * 2.0 / 100.0;
            std::complex<double > a = std::complex<double > (x, y);
            w1 = wofz(a);
            //w2 = w(a);
            //BOOST_CHECK_CLOSE(w1.real(), w2.real(), tolerance);
            //BOOST_CHECK_CLOSE(w1.imag(), w2.imag(), tolerance);
        }
    }
}
