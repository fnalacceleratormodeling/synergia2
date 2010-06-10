////////////////////////////////////////////////////////////
//
// File:          code_3_1_2.cc
// Author:        Leo Michelotti
//
// Revision date: December 6, 2007
//
////////////////////////////////////////////////////////////
//
// * Builds ERTML model from ertml_filecalls.xsif
//
// * Initial bunch is generated randomly rather than
//   read from Amundson's file, ertml_beam.dat.
//   : struct Parameters introduced to allow, in principle,
//     changing emittances on command line, though in fact
//     THIS IS NOT IMPLEMENTED.
//   : emittances taken from GDR document on RTML.
//   : bunch distribution calculated according to alphas and
//     betas given in egetaway.xsif.
//   : eigenvector (i.e. linear normal form) method used
//     : eigenvector matrix set up by hand
//       : transverse pieces from alphas and betas
//       : longitudinal piece in simplest manner
//   : uniform distribution used for action coordinates
//     : thus, distribution is NOT Gaussian. (!)
//
// * Propagates bunch element by element
//   : output ( s, "emittance" ) after each element.
//
// * Introduction of (*RIDICULOUS*) function calcTheNumber
//   to compute "emittance" = sqrt det state covariance
//
// * Four "emittances" are calculated and printed
//   : 6x6, 4x4, 2x2H, and 2x2V.
//   : calcTheNumber modified in *RIDICULOUS* fashion to handle this.
//   : output is sqrt det covariance divided by its initial value.
//     This bypasses the issue of normalization.
//
////////////////////////////////////////////////////////////

#include <parsers/xsif/XSIFFactory.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <iomanip>

#include <beamline/Particle.h>
#include <beamline/JetParticle.h>
#include <beamline/ParticleBunch.h>
#include <beamline/RefRegVisitor.h>

#include <fstream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

typedef boost::uniform_real<double>
basUnifGen;
typedef boost::variate_generator<boost::minstd_rand&, boost::uniform_real<double> >
varUnifGen;
typedef boost::normal_distribution<double>
basGaussGen;
typedef boost::variate_generator<boost::minstd_rand&, boost::normal_distribution<double> >
varGaussGen;


extern beamline* DriftsToSlots( beamline const& argbml );


enum CalcMode { m6x6 = 0, m4x4, m2x2h, m2x2v };
enum { i_x = 0, i_y, i_t, i_px, i_py, i_dpp };


struct Parameters
{
    // All parameters represent initial conditions
    double eps_1;    // [pi mm-mr]  max horizontal invariant emittance
    double eps_2;    // [pi mm-mr]  max vertical   invariant emittance
    double eps_3;    // [eV-sec]    max longitudinal emittance
    double dx;       // [mm]        horizontal beam offset
    double dy;       // [mm]        vertical   beam offset
    double dt;       // [musec]     injection time error: dt > 0 -> late
    int    n;        //             number of protons in the bunch
};


bool getParameters( Parameters& w, int argc, char** argv )
{
    w.eps_1 = 8.0e-6;
    w.eps_2 = 20.0e-9;
    w.eps_3 = 1.35e-4;
    w.dx    = 0;
    w.dy    = 0;
    w.dt    = 0;
    w.n     = 100000;

    return true;
}


double calcTheNumber( ParticleBunch& bunch, CalcMode mode )
{
    static Matrix cov(6, 6);
    static Matrix const zeroMatrix(6, 6);

    if ( m6x6 == mode ) {
        Vector centroid(6);
        double n = bunch.size();
        for ( ParticleBunch::iterator it = bunch.begin();
                it != bunch.end();
                ++it                ) {
            centroid += (*it).State();
        }
        centroid /= n;

        cov = zeroMatrix;
        for ( ParticleBunch::iterator it = bunch.begin();
                it != bunch.end();
                ++it                ) {
            for ( int i = 0; i < 6; i++ ) {
                for ( int j = 0; j < 6; j++ ) {
                    cov(i, j) +=   ( (*it).State()[i] - centroid[i] )
                                   * ( (*it).State()[j] - centroid[j] );
                }
            }
        }
        cov = cov / n;

        return sqrt(cov.determinant());
    } else if ( m4x4 == mode ) {
        int other [] = { 0, 3, 1, 4 };
        Matrix M(4, 4);
        for ( int i = 0; i < 4; i++ ) {
            for ( int j = 0; j < 4; j++ ) {
                M(i, j) = cov( other[i], other[j] );
            }
        }

        return sqrt(M.determinant());
    } else if ( m2x2h == mode ) {
        int other [] = { 0, 3 };
        Matrix M(2, 2);
        for ( int i = 0; i < 2; i++ ) {
            for ( int j = 0; j < 2; j++ ) {
                M(i, j) = cov( other[i], other[j] );
            }
        }

        return sqrt(M.determinant());
    } else if ( m2x2v == mode ) {
        int other [] = { 1, 4 };
        Matrix M(2, 2);
        for ( int i = 0; i < 2; i++ ) {
            for ( int j = 0; j < 2; j++ ) {
                M(i, j) = cov( other[i], other[j] );
            }
        }

        return sqrt(M.determinant());
    }

    return 0.0;
}

void
save_em(ParticleBunch &pb, int save_count)
{
    std::stringstream complicated;
    complicated << "code_parts_" << save_count << ".dat";
    ofstream out(complicated.str().c_str());
    int count = 0;
    for (ParticleBunch::iterator it = pb.begin(); count < 10; ++it) {
        if (count == 0) {
            std::cout << &(*it) << std::endl;
        }
        ++count;
        Vector s = it->State();
        out << std::setprecision(20);
        out << s[0] << " ";
        out << s[3] << " ";
        out << s[1] << " ";
        out << s[4] << " ";
        out << s[2] << " ";
        out << s[5] << std::endl;
    }
    out.close();
}

int main( int argc, char** argv )
{
    // Read the parameters from the command line
    // -----------------------------------------
    Parameters param;
    if ( !getParameters( param, argc, argv ) ) {
        cout << "WHOOPS!" << endl;
        return 4;
    }

    std::string file_name("ertml_filecalls.xsif");
    std::string lattice_name("ERTML");
    XSIFFactory xfactory(file_name);
    // OLD: BmlPtr bmln( DriftsToSlots( *(xfactory.create_beamline( lattice_name, xfactory.getBrho() )) ) );
    BmlPtr bmln( DriftsToSlots( *(xfactory.create_beamline( lattice_name )) ) );

    double const energy = bmln->Energy();
    Positron positron( energy );
    RefRegVisitor rrv( positron );
    bmln->accept( rrv );

    //
    // --------------------------------------
    // !XX  TWISS_EGETAWAY : BETA0, BETX = 1.6479, ALFX =  0.4982, &
    // !XX                          BETY = 8.8630, ALFY = -2.3771

    double const BETX =  1.6479;
    double const BETY =  8.8630;
    double const ALFX =  0.4982;
    double const ALFY = -2.3771;

    MatrixC E(6, 6);
    E( i_x,  i_x ) = sqrt( BETX / 2.0 );
    E( i_px, i_x ) = std::complex<double>( - ALFX / sqrt( 2.0 * BETX ), - 1.0 / sqrt( 2.0 * BETX ) );
    E( i_y,  i_y ) = sqrt( BETY / 2.0 );
    E( i_py, i_y ) = std::complex<double>( - ALFY / sqrt( 2.0 * BETY ), - 1.0 / sqrt( 2.0 * BETY ) );

    E( i_x,  i_px ) = sqrt( BETX / 2.0 );
    E( i_px, i_px ) = std::complex<double>( - ALFX / sqrt( 2.0 * BETX ),  1.0 / sqrt( 2.0 * BETX ) );
    E( i_y,  i_py ) = sqrt( BETY / 2.0 );
    E( i_py, i_py ) = std::complex<double>( - ALFY / sqrt( 2.0 * BETY ),  1.0 / sqrt( 2.0 * BETY ) );

    E( i_t,   i_t   ) = std::complex<double>( 1.0 / sqrt(2.0),    0.0           );
    E( i_dpp, i_t   ) = std::complex<double>( 0.0,            - 1.0 / sqrt(2.0) );
    E( i_t,   i_dpp ) = std::complex<double>( 1.0 / sqrt(2.0),    0.0           );
    E( i_dpp, i_dpp ) = std::complex<double>( 0.0,              1.0 / sqrt(2.0) );
    // --------------------------------------
    //

    //
    // --------------------------------------
    Positron samplePositron( energy );
    MatrixC a(6, 1);
    double psi;
    ParticleBunch bunch(samplePositron);

    Vector state(6);

    boost::minstd_rand generator(42u);
    basUnifGen angleRan(-M_PI, M_PI);
    basUnifGen xActionRan(0.0, param.eps_1 / (M_TWOPI));
    basUnifGen yActionRan(0.0, param.eps_2 / (M_TWOPI));
    basUnifGen zActionRan(0.0, param.eps_3 / (M_TWOPI));
    varUnifGen psiDist(generator, angleRan);
    varUnifGen xActionDist(generator, xActionRan);
    varUnifGen yActionDist(generator, yActionRan);
    varUnifGen zActionDist(generator, zActionRan);

    for ( int i = 0; i < 1000; i++ ) {
        psi = psiDist();
        a(0) = sqrt( xActionDist() ) * std::complex<double>( -sin(psi), cos(psi) );
        a(3) = std::conj(a(0));

        psi = psiDist();
        a(1) = sqrt( yActionDist() ) * std::complex<double>( -sin(psi), cos(psi) );
        a(4) = std::conj(a(1));

        psi = psiDist();
        a(2) = sqrt( zActionDist() ) * std::complex<double>( -sin(psi), cos(psi) );
        a(5) = std::conj(a(2));

        state = real( E * a );
        samplePositron.State() = state;
        bunch.append( samplePositron );
    }

    double const o6x6  = calcTheNumber( bunch, m6x6 );
    double const o4x4  = calcTheNumber( bunch, m4x4 );
    double const o2x2h = calcTheNumber( bunch, m2x2h );
    double const o2x2v = calcTheNumber( bunch, m2x2v );
    
    int order = 2;
    JetParticle::createStandardEnvironments(order);
    save_em(bunch, -1);
    cout << "0.0  1.0  1.0  1.0  1.0" << endl;
    double s = 0.0;
    int step = 0;
    for ( beamline::deep_iterator it = bmln->deep_begin();
            it != bmln->deep_end();
            ++it                       ) {
        JetPositron jet( energy );
        (*it)->propagate(jet);
        Mapping map(jet.State());
        for ( ParticleBunch::iterator pbit = bunch.begin();
                pbit != bunch.end();
                ++pbit) {
            Vector state(pbit->State());
            state = map(state);
            pbit->State() = state;
        }
        save_em(bunch, step);
        s += (*it)->OrbitLength(positron);
        //~ cout << s << "  " << calcTheNumber( bunch, m6x6 )/o6x6
        //~ << "  " << calcTheNumber( bunch, m4x4 )/o4x4
        //~ << "  " << calcTheNumber( bunch, m2x2h)/o2x2h
        //~ << "  " << calcTheNumber( bunch, m2x2v)/o2x2v
        //~ << endl;
        cout << step << " " << (*it)->Name() << endl;
        ++step;
    }
    // --------------------------------------
    //

    return 0;
}
