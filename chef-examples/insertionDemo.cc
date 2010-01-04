////////////////////////////////////////////////////////////
// 
// File:          insertionDemo.cc
// Author:        Leo Michelotti
// 
// Original date: March 26, 2007.  (Monday)
// Revision date: March 11, 2008.  (Tuesday)
// 
// Demonstration program written for Eric Stern.
// 
// This program has horrible memory leaks, but I don't care;
// I really don't.
// 
////////////////////////////////////////////////////////////

#include <beamline/beamline.h>
#include <beamline/marker.h>
#include <beamline/Particle.h>
#include <bmlfactory/MAD8Factory.h>


int main( int argc, char** argv ) 
{
  // -----------------------------
  // Process the command line ...
  // -----------------------------
  if( argc < 4 ) {
    cout << "Usage: " << argv[0] << " <filename>  <linename>  <partition size>" << endl;
    return -1;
  }
  cout << "Command line: ";
  for( int i = 0; i < argc; i++ ) {
    cout << argv[i] << "  ";
  }
  cout << '\n' << endl;

  std::string filename( argv[1] );
  std::string linename( argv[2] );
  int M = atoi( argv[3] );


  // --------------------------------
  // Have the factory create a model
  // --------------------------------
  MAD8Factory factory( filename );
  // REMOVE: BmlPtr bmlPtr = factory.create_beamline( linename, factory.getBrho() );
  BmlPtr bmlPtr( factory.create_beamline( linename ) );
  cout << "Energy: " << bmlPtr->Energy() << endl;
  cout << "Energy: " << factory.getEnergy() << endl;
  double energy = 1.33828;
  bmlPtr->setEnergy(energy);


  // ------------------------------------------
  // Break the model into equally spaced pieces
  // ------------------------------------------
  Proton pr(energy);
  double markerSpacing = bmlPtr->OrbitLength(pr) / ((double) M);
  std::list<std::pair<ElmPtr,double> > theMarkers;
  MarkerPtr mmm( new marker("the mark"));

  double s = 0;
  for( int i = 0; i < M-1; i++ ) {
    s += markerSpacing;
    theMarkers.push_back(std::make_pair(mmm,s));
  }

  s = 0;
  bmlPtr->InsertElementsFromList( pr, s, theMarkers );
  bmlPtr->append( mmm );


  // -----------------------------
  // Print elements in the line
  // -----------------------------
  s = 0;
  for ( beamline::deep_iterator it = bmlPtr->deep_begin(); 
        it != bmlPtr->deep_end(); 
        ++it ) {
    s += (*it)->Length();
    cout << s << ": " 
         << (*it)->Type() << "  " 
         << (*it)->Name() 
         << endl;
  }


  // -----------------------------
  // Exit
  // -----------------------------

  // I'm not bothering with garbage collection
  // and memory leaks in this program.

  return 0;
}
