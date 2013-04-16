#ifndef COMPLEX_ERROR_FUNCTION_H_
#define COMPLEX_ERROR_FUNCTION_H_

#include <cmath>
#include <complex>
#include "synergia/foundation/math_constants.h"

///////////////////////////////////////////////////////////////////////////////
//
// csp: This code is lifted from toms/680 subroutine WOFZ in www.netlib.org
// Given a complex number z, this function computes the value of Faddeeva 
// function w(z) = exp(-z^2) * erfc (- i * z), where efc is the complex 
// complementary error function and i means sqrt(-1).
// The accuracy of the algorithm for z in the 1st and 2nd quadrant is 14
// significant digits; in the 3rd and 4th, it is 13 significant digits outside
// a circular region with radius 0.126 around a zero of the function.
//
///////////////////////////////////////////////////////////////////////////////

typedef std::complex<double > Complex;

inline int
nearest_int_floor(const double x)
{
    int ix = x >= 0. ? static_cast<int > (floor(x + 0.5)) : static_cast<int > (-floor(0.5 - x));
    return ix;
}

inline Complex
wofz(Complex z)
{
    const double factor = 2.0 / sqrt(mconstants::pi);
    const double rmaxreal(0.5e154);
    const double rmaxexp(708.503061461606);
    const double rmaxgoni(3.53711887601422e15);

    double xabs = std::abs(z.real());
    double yabs = std::abs(z.imag());
    double x = xabs / 6.3;
    double y = yabs / 4.4;

    if ((xabs > rmaxreal) || (yabs > rmaxreal)) {
        throw std::runtime_error("overflow of (x^2 + y^2) in complex_error_function");
    }

    double qrho = x * x + y * y;
    double xquad = xabs * xabs - yabs * yabs;
    double yquad = 2.0 * xabs * yabs;

    double u=0.0, v=0.0;
    double u1, u2 = 0.0, v1, v2=0.0;

    if (qrho < 0.085264) {
    // If (qrho < 0.085264) then the Faddeeva function is evaluated using a 
    // power series (Abaramowitz and Stegun, Equation (7.1.5), p297)
    // n is the minimum number of terms needed to obtain the required accuracy.
        qrho = (1.0 - 0.85 * y) * sqrt(qrho);
        int n = nearest_int_floor(6.0 + 72.0 * qrho);
        int j = 2 * n + 1;
        double xsum = 1.0 / j;
        double ysum = 0.0;
        for (int i = n; i > 0; --i) {
            j -= 2;
            double xaux = (xsum * xquad - ysum * yquad) / i;
            ysum = (xsum * yquad + ysum * xquad) / i;
            xsum = xaux + 1.0 / j;
        u1 = - factor * (xsum * yabs + ysum * xabs) + 1.0;
        v1 = factor * (xsum * xabs - ysum * yabs);
        double daux = std::exp(-xquad);
        u2 = daux * cos(yquad);
        v2 = -daux * sin(yquad);

        u = u1 * u2 - v1 * v2;
        v = u1 * v2 + v1 * u2;
        }
    } else {
        double h, h2=0.0, qlambda=0.0;
        int kapn, nu;
        if (qrho > 1.0) {
        // If (qrho > 1.0) then w(z) is evaluated using the Laplace Continued 
        // Fraction
        // nu is the minimum number of terms needed to obtatin the required 
        // accuracy
            h = 0.0;
            kapn = 0;
            qrho = sqrt(qrho);
            nu = nearest_int_floor(3.0 + (1442.0 / (26.0 * qrho + 77.0)));
        } else {
        // If (0.085264 < qrho < 1.0) then w(z) is evaluated by a truncated 
        // Taylor expansion, where the Lapace Continued Fraction is used to 
        // calculate the derivatives of w(z).
        // kapn is the minimum number of terms in the Taylor expansion needed
        // to obtain the required accuracy.
        // nu is the minimum number of terms of the Continued Fraction needed 
        // to calculated the derivatives with the required accuracy.
            qrho = (1.0 - y) * sqrt(1.0 - qrho);
            h = 1.88 * qrho;
            h2 = 2.0 * h;
            kapn = nearest_int_floor(7.0 + 34.0 * qrho);
            nu = nearest_int_floor(16.0 + 26.0 * qrho);
        }

        if (h > 0.0) qlambda = std::pow(h2, kapn);

        double rx = 0;
        double ry = 0;
        double sx = 0;
        double sy = 0;

        for (int i = nu; i >= 0; --i) {
            int np1 = i + 1;
            double tx = yabs + h + np1 * rx;
            double ty = xabs - np1 * ry;
            double c = 0.5 / (tx * tx + ty * ty);
            rx = c * tx;
            ry = c * ty;
            if ((h > 0.0) && (i <= kapn)) {
                tx = qlambda + sx;
                sx = rx * tx - ry * sy;
                sy = ry * tx + rx * sy;
                qlambda = qlambda / h2;
            }
        }
        if (h == 0.0) {
            u = factor * rx;
            v = factor * ry;
        } else {
            u = factor * sx;
            v = factor * sy;
        }
        if (yabs == 0.0) u = std::exp(-xabs * xabs);
    }

    if (z.imag() < 0.0) {
        if (qrho < 0.085264) {
            u2 = 2.0 * u2;
            v2 = 2.0 * v2;
        } else {
            xquad = -xquad;

            if ((yquad > rmaxgoni) || (xquad > rmaxexp)) {
                throw std::runtime_error("overflow of (2.*exp(-z^2)) in complex_error_function");
            }

            double w1 = 2.0 * std::exp(xquad);
            u2 = w1 * cos(yquad);
            v2 = - w1 * sin(yquad);
        }

        u = u2 - u;
        v = v2 - v;
        if (z.real() > 0.0) v = -v;
    } else {
        if (z.real() < 0.0) v = -v;
    }

    return Complex (u, v);
}

//////////////////////////////////////////////////////////////
//
// Leo's version
//
// See basic_toolkit/src/erf.cc
//
//////////////////////////////////////////////////////////////

Complex wofz2( Complex );

//////////////////////////////////////////////////////////////

Complex erf( const Complex& z )
{   
    if ( ( std::abs(imag(z)) > 3.9 ) || ( std::abs(real(z)) > 3.0 ) ) {
        Complex u( - imag(z), real(z) );
        return ( 1.0 - std::exp(u*u)*wofz2(u) );
    }

    static Complex series;
    static Complex oldseries;
    static Complex arg;
    static Complex term;
    static double  den;
    static double  fctr_x;

    series        = 1.0;
    oldseries     = 0.0;
    arg           = - z * z;
    den           = 1.0;
    term          = 1.0;
    fctr_x        = 0.0;

    while ( series != oldseries ) {
        oldseries = series;
        den      += 2.0;
        fctr_x   += 1.0;
        term     *= arg / fctr_x;
        series   += term / den;
    }

    return (2.0 / sqrt(mconstants::pi)) * z * series;
}

//////////////////////////////////////////////////////////////

Complex wofz2( Complex z )
{
    static const Complex mi( 0., -1. );
    static double x;
    static double y;

    x = real(z);
    y = imag(z);

    if ( y < 0.0 )
        return 2.0*std::exp( -z*z ) - wofz2( -z );

    if ( x < 0.0 )
        return conj( wofz2( Complex( - x, y ) ) );

    if ( ( x > 6.0 ) || ( y > 6.0 ) )
        return ( - mi * z * (
                     ( 0.5124242  / ( z*z - 0.2752551 )) +
                     ( 0.05176536 / ( z*z - 2.724745  ))
                 )
               );

    if ( ( x > 3.9 ) || ( y > 3.0 ) )
        return ( - mi * z * (
                     ( 0.4613135   / ( z*z - 0.1901635 )) +
                     ( 0.09999216  / ( z*z - 1.7844927 )) +
                     ( 0.002883894 / ( z*z - 5.5253437 ))
                 )
               );

    return std::exp( -z*z )*( 1.0 - erf( mi*z ) );

}

#endif /* COMPLEX_ERROR_FUNCTION_H_ */
