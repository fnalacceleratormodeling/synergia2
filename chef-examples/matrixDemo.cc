////////////////////////////////////////////////////////////
// 
// File:          matrixDemo.cc
// Author:        Leo Michelotti
// Original date: March 20, 2007.  (Tuesday)
// Revision date: March 11, 2008.  (Tuesday)
// 
// This program has horrible memory leaks, but I don't care.
// 
////////////////////////////////////////////////////////////
// 
// Demonstration program written for Eric Stern.
// 
// Constructs a simple beamline model, breaks it into
// equally spaced pieces, and constructs transfer matrices
// for each piece.
// 
// Command line: matrixDemo <options>  N  M
// where N = number of FODO cells in the model;
//       M = number of pieces the model will be broken into.
// 
// Instructions for options are obtained by invoking
// the program with no arguments.
// 
////////////////////////////////////////////////////////////

#include <mxyzptlk/Mapping.h>
#include <beamline/quadrupole.h>
#include <beamline/drift.h>
#include <beamline/marker.h>
#include <beamline/beamline.h>
#include <beamline/Particle.h>
#include <beamline/JetParticle.h>

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

struct Options
{
  double length_F;
  double length_D;
  double abs_focalLength_F;
  double abs_focalLength_D;
  double quadCenterSpacing;
  double totalEnergy;

  static const char* usageMessage;

  Options( int, char**, int );
};


const char* Options::usageMessage = 
"  <options>  N  M"
"\nwhere N = number of FODO cells in the model."
"\n      M = number of pieces the model will be broken into."
"\n"
"\nOptions are as follows:"
"\n-------------------------------------------------------"
"\nFlag  Units   Interpretation                   Defaults"
"\n-------------------------------------------------------"
"\n-lf   meters  length of F quad              	  0.2"
"\n-ld   meters  length of D quad              	  0.2"
"\n-flf  meters  focal length of F quad        	  7"
"\n-fld  meters  focal length of D quad        	  7"
"\n-d    meters  distance between quad centers 	 10"
"\n-e    GeV     total energy of proton        	500"
"\n-ke   GeV     kinetic energy of proton"
"\n-p    GeV/c   momentum of proton"
"\n"
"\nFor example: \"-lf 0.1 -p 800\" indicates a F quad focal length"
"\nof 10 cm and proton momentum of 800 GeV/c.";


Options::Options( int argc, char** argv, int lastargs )
: length_F(0.2)
, length_D(0.2)
, abs_focalLength_F(7)
, abs_focalLength_D(7)
, quadCenterSpacing(10)
, totalEnergy(500)
{
  int limit = argc-lastargs;
  string s;
  int i = 1;
  while( i < limit ) {
    s.assign( argv[i++] );
    if( '-' == s[0] ) {
      s.assign( s.substr(1) );
      if( s == "lf" ) {
        if( i < limit ) { length_F = atof( argv[i++] );
                          if( length_F < 0.05 )
                          {   length_F = 0.05; }
        }
      }
      else if( s == "ld" ) {
        if( i < limit ) { length_D = atof( argv[i++] );
                          if( length_D < 0.05 )
                          {   length_D = 0.05; }
        }
      }
      else if( s == "flf" ) {
        if( i < limit ) { abs_focalLength_F = atof( argv[i++] );
                          if( abs_focalLength_F < 3.0*length_F )
                          {   abs_focalLength_F = 3.0*length_F; }
        }
      }
      else if( s == "fld" ) { 
        if( i < limit ) { abs_focalLength_D = atof( argv[i++] );
                          if( abs_focalLength_D < 3.0*length_D )
                          {   abs_focalLength_D = 3.0*length_F; }
        }
      }
      else if( s == "d" ) { 
        if( i < limit ) { quadCenterSpacing = atof( argv[i++] );
                          if( quadCenterSpacing < 0.001 + (length_D+length_F)/2.0 )
                          {   quadCenterSpacing = 0.001 + (length_D+length_F)/2.0;  }
        }
      }
      else if( s == "e" )  { 
        if( i < limit ) { totalEnergy = atof( argv[i++] );
                          if( totalEnergy < 1.2*PH_NORM_mp ) 
                          {   totalEnergy = 1.2*PH_NORM_mp;  }
        }
      }
      else if( s == "ke" ) { 
        if( i < limit ) { totalEnergy = PH_NORM_mp + atof( argv[i++] );
                          if( totalEnergy < 1.2*PH_NORM_mp ) 
                          {   totalEnergy = 1.2*PH_NORM_mp;  }
        }
      }
      else if( s == "p" )  { 
        if( i < limit ) { double pc = atof( argv[i++] );
                          totalEnergy = sqrt(PH_NORM_mp*PH_NORM_mp + pc*pc);
                          if( totalEnergy < 1.2*PH_NORM_mp ) 
                          {   totalEnergy = 1.2*PH_NORM_mp;  }
        }
      }
      else {
        cout << "\n*** ERROR *** Unrecognized option: " << s << endl;
      }
    }
  }
}


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

int main( int argc, char** argv ) 
{
  // -----------------------------
  // Process the command line ...
  // -----------------------------
  if( argc < 3 ) {
    cout << "Usage: " << argv[0] << Options::usageMessage << endl;
    return -1;
  }
  cout << "Command line: ";
  for( int i = 0; i < argc; i++ ) {
    cout << argv[i] << "  ";
  }
  cout << '\n' << endl;


  // -----------------------------
  // Contruct the beamline model
  // -----------------------------
  int N = atoi( argv[argc-2] ); // number of FODO cells in the model
  if( N < 1 ) { N = 1; }
  int M = atoi( argv[argc-1] ); // partition number
  if( M < 2 ) { M = 2; }        // : i.e. number of equally spaced breaks

  Options myOptions( argc, argv, 2 );

  double energy    =  myOptions.totalEnergy;
  Proton pr( energy );

  double sep         = myOptions.quadCenterSpacing;
  double length_F    = myOptions.length_F;
  double length_D    = myOptions.length_D;

  // The focal lengths will be approximate.
  double fl_F        = myOptions.abs_focalLength_F;
  double fl_D        = myOptions.abs_focalLength_D;

  double grad_F      =   ( pr.ReferenceBRho() / fl_F )/length_F;
  double grad_D      = - ( pr.ReferenceBRho() / fl_D )/length_D;

  #if 0
  QuadrupolePtr F( new quadrupole( "F" , length_F, grad_F ) );
  QuadrupolePtr D( new quadrupole( "D" , length_D, grad_D ) );
  DriftPtr     O1( new drift( "O1", (sep - length_F)/2.0 ) );
  DriftPtr     O2( new drift( "O2", sep - ((length_F+length_D)/2.0) ) );
  DriftPtr     O3( new drift( "O3", (sep - length_D)/2.0 ) );
  #endif

  #if 1
  quadrupole F( "F" , length_F, grad_F );
  quadrupole D( "D" , length_D, grad_D );
  drift     O1( "O1", (sep - length_F)/2.0 );
  drift     O2( "O2", sep - ((length_F+length_D)/2.0) );
  drift     O3( "O3", (sep - length_D)/2.0 );
  #endif

  // ----------

  beamline cell;
           cell.append( O1 );
           cell.append( F  );
           cell.append( O2 );
           cell.append( D  );
           cell.append( O3 );

  beamline model;
  for( int i = 0; i < N; i++ ) {
    model.append( BmlPtr( cell.Clone() ) );
  }


  // -----------------------------
  // Break the model into pieces
  // -----------------------------
  double markerSpacing = model.OrbitLength(pr) / ((double) M);
  std::list<std::pair<ElmPtr,double> > theMarkers;
  MarkerPtr mmm( new marker("the mark"));

  double s = 0;
  for( int i = 0; i < M-1; i++ ) {
    s += markerSpacing;
    theMarkers.push_back(std::make_pair(mmm,s));
  }

  s = 0;
  model.InsertElementsFromList( pr, s, theMarkers );
  model.append( mmm );


  // -----------------------------
  // Contruct and print transfer matrices
  // -----------------------------
  JetParticle::createStandardEnvironments();
  JetProton jpr( energy );
  Mapping inState( jpr.State() );

  int step = 0;
  for ( beamline::deep_iterator it = model.deep_begin(); 
        it != model.deep_end(); 
        ++it ) {
    if( 0 == strcmp( "marker", (*it)->Type() ) ) {
      cout << "\n-------------------\nMatrix no. " << ++step << endl;
      cout << jpr.State().Jacobian() << endl;
      jpr.State() = inState;
    }
    else {
      (*it)->propagate(jpr);
    }
  }

  

  // -----------------------------
  // Exit
  // -----------------------------

  // I'm not bothering with garbage collection
  // and memory leaks in this program.

  return 0;
}
