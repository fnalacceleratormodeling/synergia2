////////////////////////////////////////////////////////////
// 
// File:          stern_demo.cc
// Author:        Leo Michelotti
// 
// REVISION HISTORY
// 
// February 16, 2009  (original version)
// 
////////////////////////////////////////////////////////////
// 
// Reads E.G.Stern's lattice file with rfcavity
// and attempts to build a sensible model from it.
// 
////////////////////////////////////////////////////////////
// 
// ------------
// COMMAND LINE
// ------------
// stern_demo [options]
// 
// ---------------
// CURRENT OPTIONS
// ---------------
// Note: N represents an integer
//       D            a  double
//       S            a  quoted string of characters
// 
// -file    S    lattice file to be read (MAD v.8 format)
//               (default: "fobodobo_s_rf.lat")
// -machine S    name of machine to be instantiated
//               (default: "model")
// -data    S    name of file to which data are written
//               (default: "stern_demo.dat")
// -phase   D    synchronous phase of the cavity  [degrees]
//               (default: 180)
// -h       N    harmonic number
//               {default: 32)
// -dpp     D    valued of dp/p used for tracking
//               (default: 0.0001
// -n       N    number of turns to track
//               (default: 256)
// 
////////////////////////////////////////////////////////////

#include <iostream>

#include <bmlfactory/MAD8Factory.h>
#include <beamline/beamline.h>
#include <beamline/Particle.h>
#include <beamline/RefRegVisitor.h>
#include <beamline/rfcavity.h>
#include <physics_toolkit/Sage.h>

// -----------------------------------
// Options ...
// -----------------------------------
struct Options
{
  std::string fileName;
  std::string machineName;
  std::string streamName;
  double      phase;          // synchronous phase [radians]
  int         h;              // harmonic number
  int         n;              // number of turns
  double      dpp;            // dp/p, used for tracking
  
  Options( int, char**, int );
};

Options::Options( int argc, char** argv, int lastargs )
:   fileName("fobodobo_s_rf.lat")
  , machineName("model")
  , streamName(std::string( argv[0] ) + std::string(".dat"))
  , phase( M_PI )
  , h(32)
  , n(256)
  , dpp(0.0001)
{
  int limit = argc-lastargs;
  string s;
  int i = 1;
  while( i < limit ) {
    s.assign( argv[i++] );
    if( '-' == s[0] ) {
      s.assign( s.substr(1) );
      if( s == "file" ) {
        if( i < limit ) {
          fileName = std::string(argv[i++]);
        }
      }
      else if( s == "machine" ) {
        if( i < limit ) {
          machineName = std::string(argv[i++]);
        }
      }
      else if( s == "phase" ) {
        if( i < limit ) {
          double xxx = atof(argv[i++]);
          while(   xxx <= -180.0 ) { xxx += 360.0; }
          while( 180.0 <     xxx ) { xxx -= 360.0; }
          phase = xxx*( M_PI/180.0 );
        }
      }
      else if( s == "h" ) {
        if( i < limit ) {
          int xxx = atoi(argv[i++]);
          if( (8< xxx) && (xxx < 1024) ) {
            h = xxx;
          }
          else {
            cerr << "*** WARNING ***  h = " << xxx << " is out of bounds"
                    "; will use " << h
                 << endl;
          }
        }
      }
      else if( s == "dpp" ) {
        if( i < limit ) {
          double xxx = atof(argv[i++]);
          if( std::abs(xxx) < 0.05 ) { 
            dpp = xxx; 
          }
          else {
            cerr << "*** WARNING ***  dp/p = " << xxx << " is out of bounds"
                    "; will use " << dpp
                 << endl;
          }
        }
      }
      else if( s == "n" ) {
        if( i < limit ) {
          int xxx = atoi(argv[i++]);
	  //          if( (0 < xxx) && (xxx < 5000) ) {
          if( 0 < xxx ) {
            n = xxx;
          }
          else {
            cerr << "*** WARNING ***  n = " << xxx << " is out of bounds"
                    "; will use " << n
                 << endl;
          }
        }
      }
      else if( s == "data" ) {
        if( i < limit ) {
          streamName = std::string(argv[i++]);
        }
      }
      else {
        cerr << "\n*** ERROR *** Unrecognized option: " << s << endl;
      }
    }
    else {
      cerr << "\n*** ERROR *** Unable to interpret command line argument: " << s << endl;
    }
  }
}



// -----------------------------
// Main program
// -----------------------------

int main( int argc, char** argv )
{
  Options myOptions( argc, argv, 0 );
  std::string progname( argv[0] );


  // Construct/instantiate a model
  // -----------------------------
  MAD8Factory factory( myOptions.fileName );

  double const brho             = factory.getBrho();
  double const central_momentum = brho*PH_CNV_brho_to_p;
  double momentum               = central_momentum;

  cout << progname << ": brho     = " << brho     << "   T-m" << endl;
  cout << progname << ": momentum = " << momentum << " Gev/c" << endl;

  BmlPtr bmlPtr = BmlPtr( new beamline );
  bmlPtr = factory.create_beamline( myOptions.machineName );

  cout << progname << ": Beamline "
       << bmlPtr->Name()               << " created with "
       << bmlPtr->countHowMany()       << " top level elements and "
       << bmlPtr->countHowManyDeeply() << " total elements."
       << endl;

  if( !(Sage::isRing( bmlPtr )) ) {
    cout << "Whoops!" << endl;
    return -1;
  }


  // Fix the cavities' attributes
  // ----------------------------
  Proton probe;
  probe.SetReferenceMomentum( momentum );
  probe.setStateToZero();  // Not really necessary here
                           // ; probe was just instantiated.

  cout << "\nFixing cavity attributes." << endl;
  int n = 0;
  for(   beamline::deep_iterator it = bmlPtr->deep_begin()
       ; it != bmlPtr->deep_end()
       ; ++it )
  {
    if( std::string("rfcavity") == std::string((*it)->Type()) ) 
    {
      boost::dynamic_pointer_cast<rfcavity>(*it)->setHarmonicNumber( myOptions.h );
      // For now, this does nothing useful. 
      // Put here as a placeholder for the future.

      boost::dynamic_pointer_cast<rfcavity>(*it)->setPhi( myOptions.phase );

      probe.SetReferenceMomentum( momentum );  // just paranoia; this already was set
      probe.setStateToZero();                  // just paranoia; this already was set
      double circumference = bmlPtr->OrbitLength( probe );
      double freq          = double(myOptions.h) * probe.Beta() * PH_MKS_c / circumference;
      boost::dynamic_pointer_cast<rfcavity>(*it)->setFrequency( freq );

      cout << "Cavity " << ++n 
           << "\n: strength  = " << (*it)->Strength()
           << "\n: length    = " << (*it)->Length()
           << "\n: phase     = " << boost::dynamic_pointer_cast<rfcavity>(*it)->getPhi()
           << "\n: omega     = " << boost::dynamic_pointer_cast<rfcavity>(*it)->getRadialFrequency()
           << "\n: f         = " << boost::dynamic_pointer_cast<rfcavity>(*it)->getRadialFrequency()/M_TWOPI
           << "\n: h         = " << boost::dynamic_pointer_cast<rfcavity>(*it)->getHarmonicNumber()
           << endl;
    }
  }


  // Registration
  // ------------
  probe.setStateToZero();
  RefRegVisitor registrar( probe );
  bmlPtr->accept( registrar );


  // Generate some data
  // ------------------
  ofstream dataStream( myOptions.streamName.c_str() );


  probe.setStateToZero();
  probe.set_ndp( myOptions.dpp );

  dataStream << probe.get_cdt() << "  " << probe.get_ndp() << endl;
  for( int i = 0; i < myOptions.n; ++i ) {
    bmlPtr->propagate( probe );
    dataStream << probe.get_x() << " " << probe.get_npx() << " " << probe.get_y() << " " <<
      probe.get_npy() <<" " << probe.get_cdt() << "  " << probe.get_ndp() << endl;
  }
  dataStream.close();

  return 0;
}
