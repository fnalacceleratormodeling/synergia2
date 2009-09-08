////////////////////////////////////////////////////////////
//
// File:          ext_sim_e.cc
// Author:        Leo Michelotti
//
// REVISION HISTORY
//
// January, 2009  (original version)
// * Still a mess
//   ; a great many things could be improved.
//   ; most are probably not worth the effort.
//
////////////////////////////////////////////////////////////
//
// Simulates 3rd integer extraction from the debuncher.
//
////////////////////////////////////////////////////////////
//
// ------------
// COMMAND LINE
// ------------
// ext_sim_e [options]
//
// ---------------
// CURRENT OPTIONS
// ---------------
// Note: N represents an integer
//       D            a  double
//       S            a  quoted string of characters
//
// -file    S    lattice file to be read (MAD v.8 format)
//               (default: "Debunch_modified.lat")
// -machine S    name of machine to be instantiated
//               (default: "debunch")
// -tune    D    starting tune for the resonance squeeze
//               (default: 9.63)
// -phase   D    complex phase of the resonance coupling constant
//               (default: pi/2)
// -eps     D    initial, normalized invariant emittance of bunch in mm-mr
//               (default: 20; i.e. emittance = 20 pi / (beta gamma) mm-mr)
// -tol     D    tolerance for tune control
//               (default: 1.0e-8)
// -ext     D    fraction of 16 msec used to extract
//               (default: 7.0/8.0)
// -kick    D    kick imparted by electrostatic septum (radians)
//               (default: 0.001)
// -pop     N    number of protons in the bunch
//               (default: 64)
// -septum       flag to rearrange observation point to
//               upstream of septum
// -lambertson   flag to rearrange observation point to
//               upstream of lambertson
// -flip    N    changes sign of sextupole circuits
//               ; if present, N must be 1 or 2.
//
////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <basic_toolkit/PhysicsConstants.h>
#include <bmlfactory/MAD8Factory.h>
#include <beamline/beamline.h>
#include <physics_toolkit/Sage.h>
#include <beamline/marker.h>
#include <beamline/Particle.h>
#include <beamline/JetParticle.h>
#include <beamline/RefRegVisitor.h>
#include <beamline/ParticleBunch.h>
#include <beamline/septum.h>
#include <beamline/lambertson.h>
#include <physics_toolkit/BeamlineContext.h>
#include <mxyzptlk/mxyzptlk.h>
#include <physics_toolkit/TuneAdjuster.h>
#include <physics_toolkit/BmlUtil.h>

#include "s2_fish/macro_bunch_store.h"

// ----------------------
// Modelled on files:
// ~/projects/CURRENT/ILC/software/demos/src/demo_3.cc
// ~/projects/CURRENT/people/amundson_spentzouris/october_2007/emittance_calculation/code_6.cc
// ----------------------
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

int convert_chef_index(int impact_index)
{
  return impact_index/2+3*(impact_index%2);
}

void
chef_pbunch_to_mbunch_store(ParticleBunch &chef_bunch,
    Macro_bunch_store &mbs)
{
    Vector chef_state(6);
    int partnum = 0;
    for (ParticleBunch::iterator it = chef_bunch.begin(); 
        it != chef_bunch.end(); ++it )  {
        chef_state = it->State();
        for (int impact_index=0; impact_index<6; ++impact_index) {
            int chef_index = convert_chef_index(impact_index);
            mbs.local_particles(impact_index,partnum) = chef_state[chef_index]*
                mbs.units(impact_index);
        }
        mbs.local_particles(6,partnum) = partnum;
        ++partnum;
    }
 }
    
typedef boost::uniform_real<double>
basUnifGen;
typedef boost::variate_generator<boost::minstd_rand&, boost::uniform_real<double> >
varUnifGen;
typedef boost::normal_distribution<double>
basGaussGen;
typedef boost::variate_generator<boost::minstd_rand&, boost::normal_distribution<double> >
varGaussGen;

using namespace std;

extern beamline* DriftsToSlots( beamline const& );


// -----------------------------------
// Options ...
// -----------------------------------
struct Options
{
    std::string fileName;
    std::string machineName;
    std::string ostreamName;
    double startingTune_h;
    double startingTune_v;
    double endingTune;
    double deltaTune;
    double tuneTolerance;
    double adjusterTuneStep;
    double resonantTune;
    double extractionFraction;
    double septumStrength;
    double distanceToWire;
    double width;
    double gap;
    double dppMax;
    double eps_1;    // [pi mm-mr]  max horizontal invariant emittance
    double eps_2;    // [pi mm-mr]  max vertical   invariant emittance
    double eps_3;    // [eV-sec]    max longitudinal emittance (NOT USED)
    double safetyFactor;
    int    n;        //             number of protons in the bunch
    bool   startAtSeptum;
    bool   startAtLambertson;
    bool   flip[2];
    double phase_g;

    Options( int, char**, int );
};

Options::Options( int argc, char** argv, int lastargs )
        :   fileName("Debunch_modified.lat")
        , machineName("debunch")
        , startingTune_h(9.63)
        , startingTune_v(9.75)
        , resonantTune(29.0 / 3.0)
        , tuneTolerance(1.0e-8)
        , adjusterTuneStep(0.005)
        , extractionFraction(14.0 / 16.0)
        , septumStrength(0.001)
        , distanceToWire(0.016)
        , width(0.0001)
        , gap(0.1)
        , eps_1(20.0*M_PI)
        , eps_2(20.0*M_PI)
        , eps_3( 0.0)
        , safetyFactor(1.05)
        , n(64000)
        , startAtSeptum(false)
        , startAtLambertson(false)
        , phase_g(M_PI / 2.0)
{
    ostreamName = std::string( argv[0] ) + std::string("_");
    flip[0] = false;
    flip[1] = false;

    int limit = argc - lastargs;
    string s;
    int i = 1;
    while ( i < limit ) {
        s.assign( argv[i++] );
        if ( '-' == s[0] ) {
            s.assign( s.substr(1) );
            if ( s == "file" ) {
                if ( i < limit ) {
                    fileName = std::string(argv[i++]);
                }
            } else if ( s == "machine" ) {
                if ( i < limit ) {
                    machineName = std::string(argv[i++]);
                }
            } else if ( s == "tune" ) {
                if ( i < limit ) {
                    double x = atof(argv[i++]);
                    if ( (9.62 <= x) && (x <= 9.72) ) {
                        startingTune_h = x;
                    }
                }
            } else if ( s == "phase" ) {
                if ( i < limit ) {
                    double x = atof(argv[i++]);
                    while ( x <= -180. ) { x += 360.0; }
                    while ( 180. < x   ) { x -= 360.0; }
                    phase_g = x * M_PI / 180.0;
                }
            } else if ( s == "eps" ) {
                if ( i < limit ) {
                    eps_1 = M_PI * std::abs( atof(argv[i++]) );
                }
            } else if ( s == "tol" ) {
                if ( i < limit ) {
                    tuneTolerance = std::abs(atof(argv[i++]));
                    if ( tuneTolerance < 1.0e-11 || 0.1 < tuneTolerance ) {
                        tuneTolerance = 1.0e-8;
                    }
                }
            } else if ( s == "ext" ) {
                if ( i < limit ) {
                    extractionFraction = atof(argv[i++]);
                    if ( extractionFraction <= 0.0 || 1.0 < extractionFraction ) {
                        extractionFraction = 1.0;
                    }
                }
            } else if ( s == "kick" ) {
                if ( i < limit ) {
                    double x = atof( argv[i++] );
                    if ( 0.0 < x && x < 0.01 ) { septumStrength = x; }
                }
            } else if ( s == "pop" ) {
                if ( i < limit ) {
                    int m = atoi( argv[i++] );
                    if ( 2 <= m ) { n = m; }
                }
            } else if ( s == "septum" ) {
                startAtSeptum     = true;
                startAtLambertson = false;
            } else if ( s == "lambertson" ) {
                startAtLambertson = true;
                startAtSeptum     = false;
            } else if ( s == "flip" ) {
                if ( i < limit ) {
                    s.assign( argv[i++] );
                    if ( s == "1" )         { flip[0] = !flip[0]; } else if ( s == "2" )    { flip[1] = !flip[1]; } else if ( s == "both" ) {
                        flip[0] = !flip[0];
                        flip[1] = !flip[1];
                    } else {
                        cerr << "\n*** ERROR *** Unable to interpret parameter of command line argument: flip" << endl;
                    }
                } else {
                    cerr << "\n*** ERROR *** No command line argument given for: " << s << endl;
                }
            } else {
                cerr << "\n*** ERROR *** Unrecognized option: " << s << endl;
            }
        } else {
            cerr << "\n*** ERROR *** Unable to interpret command line argument: " << s << endl;
        }
    }
}



// -----------------------------
// Extraction criterion
// -----------------------------

bool isGone( Particle const& p )
{
    for ( int i = 0; i < 6; ++i ) {
        if ( 0 == finite( p.State()[i] ) ) {
            cout << "DGN: Will remove particle: infinite state detected: " << p.State()
            << " : with tag: " << p.getTag()
            << endl;
            return true;
        }
    }

    if ( std::string::npos != (p.getTag()).find("SEPTUM") ) {
        cout << "DGN: Will remove particle with tag: " << p.getTag()
        << ": " << p.State()
        << endl;
        return true;
    }

    return false;

    // return ( 0.03 < std::abs(p.get_x()) );
}


// -----------------
// Utility functions
// -----------------
// Lifted from graphic_tune_scan.cc and modified
// ---------------------------------------------

void adjustTune(   BeamlineContext& context
                   , double finalTune_h
                   , double finalTune_v
                   , double tolerance
                   , double adjuster_tune_step
                    , std::ofstream *outstreamptr)
{
    double nu_h = context.getHorizontalFracTune();
    double nu_v = context.getVerticalFracTune();

    finalTune_h = finalTune_h - double(int(finalTune_h));
    finalTune_v = finalTune_v - double(int(finalTune_v));

    (*outstreamptr) <<   "Initial tunes: horizontal = " << nu_h
    << "\n             : vertical   = " << nu_v
    << endl;

    double const dtune_h =   ( finalTune_h - nu_h )
                           / double( 1 + int( ( std::abs(finalTune_h - nu_h) ) / adjuster_tune_step ) );
    double const dtune_v =   ( finalTune_v - nu_v )
                           / double( 1 + int( ( std::abs(finalTune_v - nu_v) ) / adjuster_tune_step ) );
    double targetTune_h = nu_h + dtune_h;
    double targetTune_v = finalTune_v;
    int    step = 1;
    while ( (10.0*tolerance) < std::abs( finalTune_h - nu_h ) ) {
        (*outstreamptr) << "\nSTEP " << step << endl;
        while ( tolerance < std::abs( targetTune_h - nu_h ) ) {
            int errorCode = context.changeTunesBy( targetTune_h - nu_h,0);
            if ( errorCode != BeamlineContext::OKAY ) {
                if ( errorCode == BeamlineContext::NO_TUNE_ADJUSTER ) {
                    (*outstreamptr) << "\n\n*** ERROR *** No tune adjuster in the context!" << endl;
                }
                exit(-1);
            }
            errorCode = context.changeTunesBy( 0, targetTune_v - nu_v );
            if ( errorCode != BeamlineContext::OKAY ) {
                if ( errorCode == BeamlineContext::NO_TUNE_ADJUSTER ) {
                    (*outstreamptr) << "\n\n*** ERROR *** No tune adjuster in the context!" << endl;
                }
                exit(-1);
            }
            nu_h = context.getHorizontalFracTune();
            nu_v = context.getVerticalFracTune();
            (*outstreamptr) <<   "Tunes: horizontal = " << nu_h
            << " = " << targetTune_h << " + (" << (nu_h - targetTune_h) << ')'
            << "\n     : vertical   = " << nu_v
            << endl;
        }
        targetTune_h += dtune_h;
        ++step;
    }
}



std::vector<std::complex<double> >
resonanceSum(   std::vector<LattFuncSage::lattFunc> const& information
                , BmlPtr const& modelPtr
                , double brho )
{
    std::complex<double> const zero( 0.0, 0.0 );
    std::complex<double> const factor( 0.0, 1.0 / ( (6.0*sqrt(2.0))*(4.0*M_PI)*brho ) );

    // Output lattice functions at harmonic sextupoles
    double nu_x          = information[ information.size() - 1 ].psi.hor / M_TWOPI;
    double circumference = information[ information.size() - 1 ].arcLength;
    double harmonic      = double( nearestInteger(3.0 * nu_x) );
    double delta         = nu_x - harmonic / 3.0;

    std::vector<std::complex<double> >   g(2);
    g[0] = zero;
    g[1] = zero;

    int i = 0;
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("thinSextupole") ) {
            LattFuncSage::lattFunc localInfo = information[i];

            double theta  = M_TWOPI * localInfo.arcLength / circumference;
            double beta_x = localInfo.beta.hor;
            double psi_x  = localInfo.psi.hor;   // Assertion: 0 < psi_x
            double phase  = 3.0 * ( psi_x - delta * theta );

            while (  M_PI <  phase ) { phase -= M_TWOPI; }
            while ( phase <= -M_PI ) { phase += M_TWOPI; }

            double phase_2 = 3.0 * psi_x;
            while (  M_PI <  phase_2 ) { phase_2 -= M_TWOPI; }
            while ( phase_2 <= -M_PI ) { phase_2 += M_TWOPI; }

            double beta32 = beta_x * sqrt(beta_x);

            std::complex<double> increment
            =   std::complex<double>( beta32, 0.0 )
                * std::complex<double>( 2.0 * (*it)->Strength(), 0.0 )
                * std::complex<double>( cos(phase), - sin(phase) );

            if ( std::string::npos != std::string((*it)->Name()).find("DDD_50") ) {
                g[0] += increment;
            } else if ( std::string::npos != std::string((*it)->Name()).find("DDD_20") ) {
                g[1] += increment;
            } else {
                cerr << "*** ERROR *** " << __FILE__ << "," << __LINE__
                << ": " << (*it)->Name()
                << ": Incorrectly named harmonic sextupole, perhaps? "
                << endl;
            }
        }
        ++i;
    }

    g[0] *= factor;
    g[1] *= factor;

    return g;
}



class Debuncher
{
private:
    int notused();
    BmlPtr modelPtr;
    ParticleBunch *particle_bunch_ptr;
    Options *options_ptr;
    ofstream *outstreamptr;
    std::string progname;
    double brho;
    char **argv;
    double momentum;
    double central_momentum;
    BmlPtr bmlPtr;
public:
    //~ Debuncher(std::list<std::string> argv);
    Debuncher(std::string outfilename,double tune_h, double tune_v);
    void complete_setup();
    BmlPtr get_beamline();
    void fill_macro_bunch_store(Macro_bunch_store &mbs);
    int get_n();
};

//~ Debuncher::Debuncher(std::list<std::string> argv_list)
Debuncher::Debuncher(std::string outfilename, double tune_h, double tune_v)
{
    int argc = 1;
    argv = new char*[1];
    char myprogname[] = "prognamestr";
    argv[0] = myprogname;
    options_ptr = new Options( argc, argv, 0 );
    options_ptr->startingTune_h = tune_h;
    options_ptr->startingTune_v = tune_v;
    
    progname=argv[0] ;
    
    outstreamptr = new ofstream(outfilename.c_str());
    // Construct the model
    // -------------------
    MAD8Factory factory( options_ptr->fileName );

    brho             = factory.getBrho();
    central_momentum = brho * PH_CNV_brho_to_p;
    momentum               = central_momentum;

    (*outstreamptr) << progname << ": brho     = " << brho     << "   T-m" << endl;
    (*outstreamptr) << progname << ": momentum = " << momentum << " Gev/c" << endl;

    bmlPtr = BmlPtr( new beamline );
    bmlPtr = factory.create_beamline( options_ptr->machineName, brho );
    modelPtr = bmlPtr; // temporary; to be updated by complete_setup
    // NOTE: The second argument may no longer be needed.

    (*outstreamptr) << progname << ": Beamline "
    << bmlPtr->Name()               << " created with "
    << bmlPtr->countHowMany()       << " top level elements and "
    << bmlPtr->countHowManyDeeply() << " total elements."
    << endl;
}

void
Debuncher::complete_setup()
{
    (*outstreamptr) << progname << ": ring test (before): " << Sage::isRing( bmlPtr ) << endl;
    modelPtr = BmlPtr( DriftsToSlots( *bmlPtr ) );
    (*outstreamptr) << progname << ": ring test (after): " << Sage::isRing( modelPtr ) << endl;

    (*outstreamptr) << progname << ": Beamline "
    << modelPtr->Name()               << " created with "
    << modelPtr->countHowMany()       << " top level elements and "
    << modelPtr->countHowManyDeeply() << " total elements."
    << endl;

    Proton probe;
    probe.SetReferenceMomentum( momentum );
    probe.setStateToZero();  // Not really necessary here.

    RefRegVisitor registrar( probe );  // ??? Should put later ???
    modelPtr->accept( registrar );


#ifdef USE_SEPTUM
    // Insert the electrostatic septum
    // -------------------------------
    ThinSeptumPtr es_septum( new thinSeptum( "E_SEPTUM" ) );
    es_septum->setStrengths( options_ptr->septumStrength, 0.0 );
    es_septum->setWire( options_ptr->distanceToWire );
    es_septum->setWireWidth( options_ptr->width );
    es_septum->setGap( options_ptr->gap );

    double arcLength = 0.;
    for ( beamline::iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        arcLength += (*it)->Length();
        if ( (*it)->Name() == std::string("E_SEPTUM") ) {
            (*outstreamptr) << "Found a septum marker at s = " << arcLength << endl;
            modelPtr->putAbove( it, es_septum );
            break;
        }
    }
#endif


#ifdef USE_LAMBERTSON
    // Insert the magnetic septum (lambertson)
    // ---------------------------------------
    ThinLambPtr lambertson( new thinLamb( "LAMBERTSON" ) );

    arcLength = 0.;
    for ( beamline::iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        arcLength += (*it)->Length();
        if ( (*it)->Name() == std::string("LAMBERTSON") ) {
            (*outstreamptr) << "Found a lambertson marker at s = " << arcLength << endl;
            modelPtr->putAbove( it, lambertson );
            break;
        }
    }
#endif


    (*outstreamptr) << "-------------------------------\n\n" << endl;


    // Rearrange the model, if desired
    // -------------------------------
    if ( options_ptr->startAtSeptum ) {
        if ( 1 == modelPtr->startAt( "E_SEPTUM" ) ) {
            (*outstreamptr)  << "*** ERROR *** Permuting operation did not work." << endl;
            exit(-1);
        }
    }

    if ( options_ptr->startAtLambertson ) {
        if ( 1 == modelPtr->startAt( "LAMBERTSON" ) ) {
            (*outstreamptr)  << "*** ERROR *** Permuting operation did not work." << endl;
            exit(-1);
        }
    }



    // Establish the Jet environment
    // (necessary before continuing)
    // (failure to do this results in "mysterious" abort.
    // --------------------------------------------------
    JetParticle::createStandardEnvironments( 1 );


    // Initiate a context
    // -------------------------------------------------------
    probe.setStateToZero();
    BeamlineContext context( probe, modelPtr );


    // Set up the tune control circuits
    // -------------------------------------------------------
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( (*it)->Name().substr(0,4) == std::string("HQF1")
                || (*it)->Name().substr(0,4) == std::string("HQF2")  ) {
            (*outstreamptr) << "added " << (*it)->Name() << " to HTuneCorrector\n";
            context.addHTuneCorrector( (*it) );
        } else if ( (*it)->Name().substr(0,4) == std::string("HQD1")
                       || (*it)->Name().substr(0,4) == std::string("HQD2") ) {
            (*outstreamptr) << "added " << (*it)->Name() << " to VTuneCorrector\n";
            context.addVTuneCorrector( (*it) );
        }
    }


    // Set up the harmonic sextupole "circuits"
    // --------------------------------------
    std::vector<ElmPtr> harmonicSextupoles[2];
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("thinSextupole") ) {
            if ( std::string::npos != std::string((*it)->Name()).find("DDD_50") ) {
                (*outstreamptr)  << ": Attaching " << (*it)->Name()
                << " to first harmonic circuit."
                << endl;
                harmonicSextupoles[0].push_back( *it );
            } else if ( std::string::npos != std::string((*it)->Name()).find("DDD_20") ) {
                (*outstreamptr)  << ": Attaching " << (*it)->Name()
                << " to second harmonic circuit."
                << endl;
                harmonicSextupoles[1].push_back( *it );
            }
        }
    }


    // Calculate target emittance from invariant emittance
    // (This is the area of a triangle tangent to a circle.)
    // ---------------------------------------------------
    double beta_gamma  = central_momentum / probe.Mass();
    double emittance_x = 1.0e-6 * options_ptr->eps_1 / beta_gamma;     // meters-radians
    double emittance_y = 1.0e-6 * options_ptr->eps_2 / beta_gamma;     // meters-radians
    emittance_x *= ( 3.0 * sqrt(3.0) / M_PI );
    emittance_x *= options_ptr->safetyFactor;


    // Set initial tune.
    // Store quad settings for initial tune.
    // ----------------------------------------
    double tune_h = options_ptr->startingTune_h;
    double tune_v = options_ptr->startingTune_v;

    // context.reset(); ??? This would make adjustTune fail.  WHY ???
    //                  ??? BECAUSE the reset deleted the context's tuneadjuster.
    //                  ??? that was instantiated when the control circuits
    //                  ??? were established.
    //                  !!! FIX THIS !!!
    adjustTune( context, tune_h, tune_v, options_ptr->tuneTolerance, 
        options_ptr->adjusterTuneStep, outstreamptr );

    (*outstreamptr) << "\n\n " << progname << ": INITIAL QUADRUPOLE SETTINGS" << endl;
    std::vector<double> initialSettings;
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("quadrupole") ) {
            (*outstreamptr) << (*it)->Name() << "  " << setprecision(12) << (*it)->Strength() << endl;
            initialSettings.push_back((*it)->Strength());
        }
    }


    // ... determine the target value of | g |
    // -----------------------------------
    double delta = tune_h - options_ptr->resonantTune;
    double target_abs_g =   std::abs(delta)
                            / sqrt( 2.0 * sqrt(3.0) * emittance_x );
    double target_phase_g = options_ptr->phase_g;

    std::complex<double> target_g =
        std::complex<double>(   target_abs_g * cos(target_phase_g)
                                , target_abs_g * sin(target_phase_g) );

    (*outstreamptr) << "\nFor tune = " << tune_h
    << ", target |g| = " << target_abs_g
    << endl;


    // ... calculate current value of g
    // --------------------------------
    // context.reset();  ??? Needed ???  ??? Produces error ???
    std::vector<LattFuncSage::lattFunc> information = context.getTwissArray();

    std::vector<std::complex<double> > g = resonanceSum(   information
                                           , modelPtr
                                           , brho );

    std::complex<double> total_g = g[0] + g[1];
    double abs_g = std::abs(total_g);
    (*outstreamptr) << "\nCurrent value: g = " << total_g
    << ", |g| = " << abs_g
    << endl;


    // ... adjust harmonic sextupole strengths ...
    // ---------------------------------------
    (*outstreamptr) << "\n\n"  << ": Setting (integrated) sextupole strengths"
    << endl;


    // ...... change polarity, if desired
    // ----------------------------------
    for ( int i = 0; i < 2; ++i ) {
        if ( options_ptr->flip[i] ) {
            for (   std::vector<ElmPtr>::iterator it =  harmonicSextupoles[i].begin()
                    ; it != harmonicSextupoles[i].end()
                    ; ++it ) {
                (*it)->setStrength( - ((*it)->Strength()) );
            }
        }
    }


    // ...... scale to target value of harmonic coupling
    // -------------------------------------------------
    (*outstreamptr) << "\n\nSetting (integrated) sextupole strengths"
    << endl;
    double ratio = target_abs_g / abs_g;
    (*outstreamptr) << "\nratio = " << ratio << endl;
    for ( int i = 0; i < 2; ++i ) {
        for (   std::vector<ElmPtr>::iterator it =  harmonicSextupoles[i].begin()
                ; it != harmonicSextupoles[i].end()
                ; ++it ) {
            (*outstreamptr) << (*it)->Type() << "  " << (*it)->Name() << ": "
            << (*it)->Strength() << " -> " << ( ratio * (*it)->Strength() );
            (*it)->setStrength(   ratio * (*it)->Strength() );
            (*outstreamptr) << ":  B''l = " << 2.0*(*it)->Strength() << " T/m" << endl;
        }
        g[i] *= ratio;
    }
    (*outstreamptr) << endl;


    // ...... rotate, if desired
    // -------------------------
    double s[2];
    double denom = imag( conj(g[1]) * g[0] );
    s[0] =   imag( conj(g[1]) * target_g ) / denom;
    s[1] = - imag( conj(g[0]) * target_g ) / denom;

    (*outstreamptr) << "\n\n"  << ": Doing rotation." << endl;
    for ( int i = 0; i < 2; ++i ) {
        for (   std::vector<ElmPtr>::iterator it =  harmonicSextupoles[i].begin()
                ; it != harmonicSextupoles[i].end()
                ; ++it ) {
            (*outstreamptr) << (*it)->Type() << "  " << (*it)->Name() << ": "
            << (*it)->Strength() << " -> " << ( s[i] * (*it)->Strength() );
            (*it)->setStrength( s[i] * (*it)->Strength() );
            (*outstreamptr) << ":  B''l = " << 2.0*(*it)->Strength() << " T/m" << endl;
        }
    }
    (*outstreamptr) << endl;


    // Set resonant tune.
    // Store quad settings for resonant tune.
    // ----------------------------------------
    tune_h = options_ptr->resonantTune;

    // context.reset(); ??? This would make adjustTune fail.  WHY ???
    //                  ??? BECAUSE the reset deleted the context's tuneadjuster.
    //                  ??? that was instantiated when the control circuits
    //                  ??? were established.
    //                  !!! FIX THIS !!!
    adjustTune( context, tune_h, tune_v, options_ptr->tuneTolerance, 
        options_ptr->adjusterTuneStep, outstreamptr );



    // Store the settings
    // -------------------------------------------
    (*outstreamptr) << "\n\nFINAL QUADRUPOLE SETTINGS" << endl;
    std::vector<double> finalSettings;
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("quadrupole") ) {
            (*outstreamptr) << (*it)->Name() << "  " << setprecision(12) << (*it)->Strength() << endl;
            finalSettings.push_back((*it)->Strength());
        }
    }



    // Just a test
    // -------------------------------------------
    int i = 0;
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("quadrupole") ) {
            (*it)->setStrength( finalSettings[i++] );
        }
    }
    (*outstreamptr) << "\n\nTEST: AFTER SETTING FINAL CONDITIONS" << endl;
    context.reset();   // NOTE: destroys context's TuneAdjuster
    double nu_h = context.getHorizontalEigenTune();
    double nu_v = context.getVerticalEigenTune();
    double resFracTune = options_ptr->resonantTune - double(int(options_ptr->resonantTune));

    (*outstreamptr) <<   "Fractional tunes: horizontal = " << nu_h
    << " = resFracTune + (" << (nu_h - resFracTune) << ')'
    << "\n     : vertical   = " << nu_v
    << endl;


    // Set initial conditions, prior to ramping
    // ----------------------------------------
    i = 0;
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("quadrupole") ) {
            (*it)->setStrength( initialSettings[i++] );
        }
    }

    (*outstreamptr) << "\n\nTEST: AFTER SETTING INITIAL CONDITIONS" << endl;
    context.reset();
    nu_h = context.getHorizontalEigenTune();
    nu_v = context.getVerticalEigenTune();
    double startingFracTune_h = options_ptr->startingTune_h - 
        double(int(options_ptr->startingTune_h));
    double startingFracTune_v = options_ptr->startingTune_v -
        double(int(options_ptr->startingTune_v));
    (*outstreamptr) <<   "Tunes: horizontal = " << nu_h
    << " = startingFracTune + (" << (nu_h - startingFracTune_h) << ')'
    << "\n     : vertical   = " << nu_v
    << endl;



    // Determine ramp time
    // ----------------------
    probe.setStateToZero();
    double const c             = PH_MKS_c;
    double const circumference = modelPtr->OrbitLength( probe );
    double const v             = c * probe.ReferenceBeta();

    (*outstreamptr) << "\n\nCircumference = " <<  circumference << endl;

    (*outstreamptr) << "Infinite energy limit:" << endl;
    (*outstreamptr) << ( circumference / c )*1.0e6 << " mu-sec" << endl;
    (*outstreamptr) << ( c / circumference ) / 1000. << "    kHz" << endl;

    (*outstreamptr) << "No. of turns in 1/15 sec = "
    << ( c / circumference ) / 15.0
    << endl;

    (*outstreamptr) << "No. of turns in 1/60 sec = "
    << ( c / circumference ) / 60.0
    << endl;

    (*outstreamptr) << "\nAt 8 GeV/c momentum:" << endl;
    (*outstreamptr) << ( circumference / v )*1.0e6 << " mu-sec" << endl;
    (*outstreamptr) << ( v / circumference ) / 1000. << "    kHz" << endl;

    (*outstreamptr) << "No. of turns in 1/60 sec = "
    << ( v / circumference ) / 60.0
    << endl;

    int turns_to_extract = (int) (options_ptr->extractionFraction * ( ( v / circumference ) / 60.0 ));

    (*outstreamptr) << "Extraction will be carried out in "
    << turns_to_extract << " turns."
    << endl;


    // ----------------------
    // Create initial distribution ...
    // ----------------------
    // Modelled on files:
    // ~/projects/CURRENT/ILC/software/demos/src/demo_3.cc
    // ~/projects/CURRENT/people/amundson_spentzouris/october_2007/emittance_calculation/code_6.cc
    // ----------------------


    // ... construct properly normalized eigenvector matrix
    // ----------------------------------------------------
    JetProton jp( probe.ReferenceEnergy() );
    modelPtr->propagate(jp);

    MatrixC E(6, 6);
    Vector phase(2);
    E = ( jp.State().Jacobian() ).Matrix::eigenVectors();
    // NOTE: it should have been possible to use .eigenVectors() rather
    // than explicitly indicating Matrix::eigenVectors(). That the latter
    // is necessary is an not understood "feature" of using templates.
    BmlUtil::normalize( E, phase );


    // ... populate a bunch according to linearly invariant distribution
    // -----------------------------------------------------------------
    probe.setStateToZero();
    boost::minstd_rand generator(42u);

    basUnifGen angleRan(-M_PI, M_PI);
    basUnifGen xActionRan(0.0, emittance_x / M_TWOPI);
    basUnifGen yActionRan(0.0, emittance_y / M_TWOPI);

    varUnifGen psiDist(generator, angleRan);
    varUnifGen xActionDist(generator, xActionRan);
    varUnifGen yActionDist(generator, yActionRan);

    MatrixC a(6, 1);
    double psi;
    Proton sampleProton( probe.ReferenceEnergy() );
    //~ ParticleBunch bunch( sampleProton );
    particle_bunch_ptr = new ParticleBunch(sampleProton);

    for ( int i = 0; i < options_ptr->n; i++ ) {
        psi = psiDist();
        a(0) = sqrt( xActionDist() ) * std::complex<double>( -sin(psi), cos(psi) );
        a(3) = std::conj(a(0));

        psi = psiDist();
        a(1) = sqrt( yActionDist() ) * std::complex<double>( -sin(psi), cos(psi) );
        a(4) = std::conj(a(1));

        a(2) = std::complex<double>( 0, 0 );
        a(5) = std::conj(a(2));
        // ??? NEXT STEP: This last must be upgraded to include
        // ??? a longitudinal distribution as well. The issue
        // ??? is to make certain that BmlUtil::normalize
        // ??? handles the longitudinal sector "correctly."
        // ??? (It probably does well enough.)

        sampleProton.State() = real( E * a );
        sampleProton.State()[2] = 0.0;
        sampleProton.State()[5] = 0.0;
        particle_bunch_ptr->append( sampleProton );
    }


    // Ramp and extract ...
    // ----------------------

    // ... Initial setting of quadrupoles
    // ----------------------------------
    i = 0;
    for ( beamline::const_iterator it = modelPtr->begin();
            it != modelPtr->end();
            ++it ) {
        if ( std::string((*it)->Type()) == std::string("quadrupole") ) {
            (*it)->setStrength( initialSettings[i++] );
        }
    }


    (*outstreamptr) << "\n\nBEGIN RAMPING\n" << endl;
}

#if 0
int
Debuncher::notused()
{
    // Output data
    // -----------
    time_t encodedTime;
    time( &encodedTime );
    tm* timestamp = localtime( &encodedTime );

    ostringstream namebuffer;
    namebuffer << options_ptr->ostreamName
    << timestamp->tm_year + 1900
    << timestamp->tm_mon + 1
    << timestamp->tm_mday << '_'
    << timestamp->tm_hour << '_'
    << timestamp->tm_min  << ".dat";
    ofstream dataStream( namebuffer.str().c_str() );
    int turnNumber = 0;
    int oldSize = particle_bunch_ptr->size();
    dataStream << turnNumber << "  " << particle_bunch_ptr->size() << endl;

#if 0
    for ( ParticleBunch::iterator it = particle_bunch_ptr->begin();
            it != particle_bunch_ptr->end();
            ++it                ) {
        dataStream << (*it).get_x() << "  " << (*it).get_npx() << endl;
    }
    dataStream << "  " << endl;
#endif


    // Begin extraction
    // ----------------
    double newStrength, epsilon;
    for ( int n = 1; n <= turns_to_extract; ++n ) {
        if ( 0 == n % 100 ) {
            (*outstreamptr) << "Turn " << n << endl;
        }

        modelPtr->propagate( particle_bunch_ptr );
        turnNumber++;

        particle_bunch_ptr->remove( isGone );

        // Output data
        // -----------
        if ( particle_bunch_ptr->size() < oldSize ) {
            dataStream << turnNumber << "  " << particle_bunch_ptr->size() << endl;
            oldSize = particle_bunch_ptr->size();
        }

#if 0
        if ( 0 == n % 100 ) {
            for ( ParticleBunch::iterator it = particle_bunch_ptr->begin();
                    it != particle_bunch_ptr->end();
                    ++it                ) {
                dataStream << (*it).get_x() << "  " << (*it).get_npx() << endl;
            }
            dataStream << "  " << endl;
        }
#endif

        // Reset quadrupoles
        // -----------------
        i = 0;
        for ( beamline::const_iterator it = modelPtr->begin();
                it != modelPtr->end();
                ++it ) {
            if ( std::string((*it)->Type()) == std::string("quadrupole") ) {
                epsilon = ((double) n) / ((double) turns_to_extract);
                newStrength = epsilon * finalSettings[i] + (1.0 - epsilon) * initialSettings[i];
                (*it)->setStrength( newStrength );
                ++i;
            }
        }
    }  // End extraction


    // Output data
    // -----------
    dataStream << turnNumber << "  " << particle_bunch_ptr->size() << endl;

#if 0
    for ( ParticleBunch::iterator it = particle_bunch_ptr->begin();
            it != particle_bunch_ptr->end();
            ++it                ) {
        dataStream << (*it).get_x() << "  " << (*it).get_npx() << endl;
    }
#endif


    (*outstreamptr) << "\n\nEND RAMPING\n" << endl;

#if 0
    ofstream xyStream( "ext_sim_d_b_xy.dat" );
    ofstream xpxStream( "ext_sim_d_b_xpx.dat" );
    ofstream ypyStream( "ext_sim_d_b_ypy.dat" );
    for ( ParticleBunch::iterator it = particle_bunch_ptr->begin();
            it != particle_bunch_ptr->end();
            ++it                ) {
        xpxStream << (*it).get_x() << "  " << (*it).get_npx() << endl;
        ypyStream << (*it).get_y() << "  " << (*it).get_npy() << endl;
        xyStream  << (*it).get_x() << "  " << (*it).get_y()   << endl;
    }
    xpxStream.close();
    ypyStream.close();
    xyStream.close();
#endif

    dataStream.close();

    return 0;
}
#endif

BmlPtr
Debuncher::get_beamline()
{
    return modelPtr;
}

void
Debuncher::fill_macro_bunch_store(Macro_bunch_store &mbs)
{
    chef_pbunch_to_mbunch_store(*particle_bunch_ptr,mbs);
}

int
Debuncher::get_n()
{
    return options_ptr->n;
}

#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include "s2_fish/container_conversions.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(debuncher)
{
    //---------------------------------------------------------------------
    // std::vector<> conversions
    //---------------------------------------------------------------------
    scitbx::boost_python::container_conversions::from_python_sequence <
    std::vector<std::string>,
    scitbx::boost_python::container_conversions::variable_capacity_policy > ();

    //~ class_<Debuncher>("Debuncher",init<std::list<std::string> >() )
    class_<Debuncher>("Debuncher",init<std::string,double,double>() )
    .def("complete_setup",&Debuncher::complete_setup)
    .def("get_beamline",&Debuncher::get_beamline)
    //~ .def("zero_sextupoles",&Debuncher::zero_sextupoles)
    .def("fill_macro_bunch_store",&Debuncher::fill_macro_bunch_store)
    .def("get_n",&Debuncher::get_n)
    ;
}


