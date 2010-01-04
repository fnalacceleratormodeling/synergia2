#include <stdio.h>
#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
typedef std::complex<double> Complex;
//#include "complexAddon.h"

const Complex complex_1( 1.0, 0.0 );
const Complex complex_0( 0.0, 0.0 );
const Complex complex_i( 0.0, 1.0 );

#include "basic_toolkit/MathConstants.h"

class BasErs_field
{
private:
    double sigma[2];
public:
    BasErs_field( double* = 0
                            /* pointer to an array containing
                               sigmax and sigmay [m] */
                );
    BasErs_field( const BasErs_field& );
    ~BasErs_field();

    char useRound;         // By default = 1
    // If 1: then round beam approximation
    // used when horizontal and vertical
    // sigmas approximately equal.

    std::vector<double> NormalizedEField( double x, double y );
    /* returns the "normalized" electric field
       in the rest frame of the bunch, in inverse
       meters.  To get the field [V/m], this must
       be multiplied by Q/(2 pi epsilon_o), where
       Q is the line density of charge [C/m] (in
       rest frame). */

    void GetSigma( double* );
};

