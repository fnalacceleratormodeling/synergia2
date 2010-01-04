/*************************************************************************
**************************************************************************
**************************************************************************
******
******  BASIC TOOLKIT:  Low level utility C++ classes.
******  Version:   4.0
******
******  File:      erf.cc
******
******  Copyright (c) 1990 Universities Research Association, Inc.
******                All Rights Reserved
******
******  Author:    Leo Michelotti
******
******             Fermilab
******             P.O.Box 500
******             Mail Stop 220
******             Batavia, IL   60510
******
******             Phone: (630) 840 4956
******             Email: michelotti@fnal.gov
******
******  Usage, modification, and redistribution are subject to terms
******  of the License and the GNU General Public License, both of
******  which are supplied with this software.
******
**************************************************************************
*************************************************************************/


/*
** Complex error function.
**
** Author: Leo Michelotti
** Date: July 7, 1995
**
*/

#ifdef __VISUAL_CPP__
#include <complex>
#include <iostream>
using std::cerr;
#else
#include <complex>
#endif
#include "complexAddon.h"
#include "MathConstants.h"

Complex w( Complex );

//////////////////////////////////////////////////////////////

Complex erfSeries( const Complex& z )
{
    static Complex series;
    static Complex oldseries;
    static Complex arg;
    static double  den;
    static Complex term;

    series        = 1.0;
    oldseries     = 0.0;
    arg           = 2.0 * z * z;
    den           = 1.0;
    term          = 1.0;

    while ( series != oldseries ) {
        oldseries = series;
        den      += 2.0;
        term     *= ( arg / den );
        series   += term;
    }

    return z*series;
}

//////////////////////////////////////////////////////////////

Complex erf( const Complex& z )
{
    if ( ( fabs(imag(z)) > 3.9 ) || ( fabs(real(z)) > 3.0 ) ) {
        Complex u( - imag(z), real(z) );
        return ( 1.0 - std::exp(u*u)*w(u) );
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

    return (2.0 / MATH_SQRTPI)*z*series;
}

//////////////////////////////////////////////////////////////

Complex erfc( const Complex& z )
{
    static const Complex one( 1.0, 0.0 );
    return ( one - erf( z ) );
}

//////////////////////////////////////////////////////////////

Complex w( Complex z )
{
    static const Complex mi( 0., -1. );
    static double x;
    static double y;

    x = real(z);
    y = imag(z);

    if ( y < 0.0 )
        return 2.0*std::exp( -z*z ) - w( -z );

    if ( x < 0.0 )
        return conj( w( Complex( - x, y ) ) );

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

