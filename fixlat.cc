#include <iostream>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "basic_toolkit/PhysicsConstants.h"
#include "beamline/beamline.h"
#include "bmlfactory/bmlfactory.h"
#include "physics_toolkit/BeamlineContext.h"
#include "beamline/RefRegVisitor.h"
#include "physics_toolkit/ClosedOrbitSage.h"

using namespace std;
extern beamline* DriftsToSlots( beamline& original );

madparser* mp = 0;

int main(int argc, char **argv)
{
  if (argc!=3) {
    cout << "usage: fixlat mad_file line_name\n";
    exit(1);
  }
  string mad_file(argv[1]);
  string line_name(argv[2]);
  double kinetic_energy = 0.4;
  double mass = PH_NORM_mp;
  double energy   = kinetic_energy + mass;
  double momentum = sqrt(energy*energy - mass*mass);

  int order = 1;

  JetParticle::createStandardEnvironments(order);
  
  double brho = (fabs(momentum))/PH_CNV_brho_to_p;
  bmlfactory* bml_fact = new bmlfactory(mad_file.c_str(), brho);
  beamline* bmline_orig = bml_fact->create_beamline(line_name.c_str());
  bmline_orig->flatten();
 
  beamline* bmline = DriftsToSlots(*bmline_orig);
  Proton jfcproton(energy);
  BeamlineContext* bmln_context = new BeamlineContext(jfcproton,bmline,false);

  JetProton jpr(energy);

  ClosedOrbitSage cos(bmline);
  cos.set_verbose();
  if (!cos.findClosedOrbit(&jpr)) {
    cout << "found a closed orbit\n";
  } else {
    cout << "failed to find a closed orbit\n";
  }


  if (bmln_context->isTreatedAsRing()) {
    cout << "treated as ring\n";
  } else {
    cout << "nope, it isn't a ring\n";
  }

  bmlnElmnt* be;
  DeepBeamlineIterator deep_beamline_iterator(bmline);

  int i=0;
  while((be = deep_beamline_iterator++)) {
    //    cout << i << " " << be->Name() << " " << be->Type() << endl;
    ++i;
  }
  cout << "success!\n";
  return(0);
}
