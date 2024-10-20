#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_particles.h"

#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/bunch_simulator.h"

#include "synergia/utils/utils.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"

int macroparticles = 16;
int spectparticles = 4;
double realparticles = 285452129459.3449; // 0.1 mA in IOTA
int nturns = 100;

Lattice
get_lattice()
{
    static std::string booster_madx(R"foo(
! Simplified Booster lattice

// From JFO 2022-12-08
// The simplified booster lattice is trivial. I have no madx lattice (you could make one
// very easily) - I just use pyorbit classes directly to instantiate a basic cell;  it is then
// replicated 24 times. I took the bending magnets lengths and 
// strengths directly from the official MADX  lattice file. 

// The basic cell is 

// d1 fmag d2 dmag d3 

// d1, d2, d3 : drifts of lengths 0.6 0.5 and 3.0 m
// fmag:  focusing      bend   L = 2.889612 m 
// dmag   defocusing bend   L = 2.889612 m

// total cell length: 19.758 m
// total ring  length = 24*19.758 = 474.20 m 

// The length, focusing strengths and curvature radius of the 
// magnets are as in the booster MADX file.  

// If you entered 1 cell correctly, you should get the periodic solution:
// bx = 33.86 m ax = 0
// by = 5.39m    ay =0 
// For 24 cells, the raw tunes are nux = 7.017 and nuy = 6.675.  You will need to tweak the nominal focusing strengths a bit to avoid resonances.


//--------- Nominal Gradient Magnet Definitions  

// EGS
// The apparent cell structure is actually:

// D1, l=0.6;
// FMAGU01;
// D2, l=0.5;
// DMAGU01;
// D3, l=6.0;
// DMAGD01;
// D4, l=0.5;
// FMAGD01;
// DR, l=0.6


! Expand structure to include corrector packages:

! Normal short straight

! drift 0.6 (end of previous cell)
! drift 0.176
! corrector package (short) l=0.168
! drift 0.256
! end of short straight

! non-RF cavity long straight
!
! drift 5.581
! correction package (long) l=0.168
! drift 0.251
!

! long straight with RF cavity
!
! drift 0.21
! rfcavity drift-cavity-drift 2.35
! drift 0.12
! rfcavity drift-cavity-drift 2.35
! 4 drifts total 0.551
! corrector packages (long) 0.168
! drift 0.251

! Corrector package:
!
! HLxx, HKICKER, l=0.024
! VLxx, VKICKER, l=0.024
! Q{S|L}xx QUADRUPOLE, l=0.024 (normal)
! MULTIPOLE Q{S|L}ERR, l=0
! MONITOR, l=0.024
! QS{S|L}xx, QUADRUPOLE, l=0.024 (skew)
! MULTIPOLE QS{S|L}ERR, l=0
! SEXTUPOLE SX{L|S}, l=0.024 (normal)
! SEXTUPOLE SS{L|S}, l=0.024 (skew)
! 

ke1 = 0.8;  !800 MeV kinetic energy at injection

rhof  :=  40.847086;   !  bending radius of focusing magnet
rhod  :=  48.034101;   !  bending radius of defocusing magnet

blength :=     2.889612;    !  arc length for both F and D magnets
blengthf :=    2.889009499; !  physical length (straight length) for F magnet
blengthd :=    2.889176299; !  physical length (straight length) for D magnet


!
! The quad field for the gradient magnet is contained in file " qsdqsf.dat" to be read in before this file !
!
! read from file at time step = 7
qsd := -57.38855012e-3;
qsf := 54.10921561e-3;


! These ssd and ssf strengths come from fitting to 01 Dec 2015 chromaticity data
! and predicts chromaticity much better than using the Drozhdin et al measurements above

ssd :=  -0.04381647074 + ke1*(0.009150934932+ ke1*(-0.0023900895  + ke1*(0.000318068028 -  ke1* 1.6353205e-05)));

ssf :=  -0.006384940088 + ke1*(0.01967542848 + ke1*( -0.006776746 + ke1*(0.00091367565 - ke1* 4.293705e-05)));

 !
 ! Gradient magnets defined by their physical length aith their bend angle
 ! being defined by the arc length/radius of curvature

!FMAG: RBEND,  L = blengthf  , ANGLE = blength/rhof, K1 = qsf  , K2 = ssf;
FMAG: SBEND,  L = blength  , ANGLE = blength/rhof, e1=blength/(2*rhof), e2=blength/(2*rhof), K1 = qsf  , K2 = ssf, type=fmag;
DMAG: SBEND,  L = blength  , ANGLE = blength/rhod, e1=blength/(2*rhod), e2=blength/(2*rhod), K1 = qsd  , K2 = ssd, type=dmag;

! drifts in the short straight section
mins: drift, l=0.5, type=shortstraight;
sc: drift, l=0.6, type=shortstraight;
sa: drift, l=0.176, type=shortstraight;
sb: drift, l=0.256, type=shortstraight;

! drifts in the long straight section
dlong: drift, l=5.581, type=longstraight; ! for no-RF cavity cells
drifta: drift, l=0.21, type=longstraight; ! start of RF cavity cell
driftb: drift, l=0.12, type=longstraight;
dmidls: drift, l=0.551, type=longstraight; ! (drift in the middle of the longstraight)
drifte: drift, l=0.251, type=longstraight; ! end of RF cavity cell

drrf: drift, l=2.35/2, type=rfaperture;
rfc: rfcavity, l=0, harmon=84, volt=0.2/22, lag=0, type=rfaperture;

! short corrector package
hsxx: drift, l=0.024, type=shortstraight;
vsxx: drift, l=0.024, type=shortstraight;
qsxx: quadrupole, k1=0.0, l=0.024, type=shortstraight;
! multipole not included
bpms: drift, l=0.024, type=shortstraight;
qssxx: drift, l=0.024, type=shortstraight; // skew quad
! multipole not included
sxsxx: sextupole, k2=0.0, l=0.024, type=shortstraight;
sssxx: drift, l=0.024, type=shortstraight; // skew sextupole

cpshort: line=(hsxx, vsxx, qsxx, bpms, qssxx, sxsxx, sssxx);

! long corrector package
hlxx: drift, l=0.024, type=longstraight;
vlxx: drift, l=0.024, type=longstraight;;
qlxx: quadrupole, k1=0.0, l=0.024, type=shortstraight;
! multipole not included
bpml: drift, l=0.024, type=shortstraight;
qlsxx: drift, l=0.024, type=longstraight; // skew quad
! multipole not included
sxlxx: sextupole, k2=0.0, l=0.024, type=longstraight;
sslxx: drift, l=0.024, type=longstraight; // skew sextupole

cplong: line=(hlxx, vlxx, qlxx, bpms, qlsxx, sxlxx, sslxx);

!!!!!!!!!!!!!!!!   beginning of ring definition

fmagu01: fmag;
fmagd01: fmag;
dmagu01: dmag;
dmagd01: dmag;

cell01 : line = (sa, cpshort, sb, fmagu01, mins, dmagu01, dlong, cplong, drifte, dmagd01, mins, fmagd01, sc);

fmagu02: fmag;
fmagd02: fmag;
dmagu02: dmag;
dmagd02: dmag;

cell02 : line = (sa, cpshort, sb, fmagu02, mins, dmagu02, dlong, cplong, drifte, dmagd02, mins, fmagd02, sc);

fmagu03: fmag;
fmagd03: fmag;
dmagu03: dmag;
dmagd03: dmag;

cell03 : line = (sa, cpshort, sb, fmagu03, mins, dmagu03, dlong, cplong, drifte, dmagd03, mins, fmagd03, sc);

fmagu04: fmag;
fmagd04: fmag;
dmagu04: dmag;
dmagd04: dmag;

cell04 : line = (sa, cpshort, sb, fmagu04, mins, dmagu04, dlong, cplong, drifte, dmagd04, mins, fmagd04, sc);

fmagu05: fmag;
fmagd05: fmag;
dmagu05: dmag;
dmagd05: dmag;

cell05 : line = (sa, cpshort, sb, fmagu05, mins, dmagu05, dlong, cplong, drifte, dmagd05, mins, fmagd05, sc);

fmagu06: fmag;
fmagd06: fmag;
dmagu06: dmag;
dmagd06: dmag;

cell06 : line = (sa, cpshort, sb, fmagu06, mins, dmagu06, dlong, cplong, drifte, dmagd06, mins, fmagd06, sc);

fmagu07: fmag;
fmagd07: fmag;
dmagu07: dmag;
dmagd07: dmag;

cell07 : line = (sa, cpshort, sb, fmagu07, mins, dmagu07, dlong, cplong, drifte, dmagd07, mins, fmagd07, sc);

fmagu08: fmag;
fmagd08: fmag;
dmagu08: dmag;
dmagd08: dmag;

cell08 : line = (sa, cpshort, sb, fmagu08, mins, dmagu08, dlong, cplong, drifte, dmagd08, mins, fmagd08, sc);
 
fmagu09: fmag;
fmagd09: fmag;
dmagu09: dmag;
dmagd09: dmag;

cell09 : line = (sa, cpshort, sb, fmagu09, mins, dmagu09, dlong, cplong, drifte, dmagd09, mins, fmagd09, sc);

fmagu10: fmag;
fmagd10: fmag;
dmagu10: dmag;
dmagd10: dmag;

cell10 : line = (sa, cpshort, sb, fmagu10, mins, dmagu10, dlong, cplong, drifte, dmagd10, mins, fmagd10, sc);

fmagu11: fmag;
fmagd11: fmag;
dmagu11: dmag;
dmagd11: dmag;

cell11 : line = (sa, cpshort, sb, fmagu11, mins, dmagu11, dlong, cplong, drifte, dmagd11, mins, fmagd11, sc);

fmagu12: fmag;
fmagd12: fmag;
dmagu12: dmag;
dmagd12: dmag;

cell12 : line = (sa, cpshort, sb, fmagu12, mins, dmagu12, dlong, cplong, drifte, dmagd12, mins, fmagd12, sc);

fmagu13: fmag;
fmagd13: fmag;
dmagu13: dmag;
dmagd13: dmag;

cell13 : line = (sa, cpshort, sb, fmagu13, mins, dmagu13, dlong, cplong, drifte, dmagd13, mins, fmagd13, sc);

fmagu14: fmag;
fmagd14: fmag;
dmagu14: dmag;
dmagd14: dmag;
rf01: line=(drrf, rfc, drrf);
rf02: line=(drrf, rfc, drrf);

cell14: line = (sa, cpshort, sb, fmagu14, mins, dmagu14, drifta, rf01, driftb, rf02, dmidls, cplong, drifte, dmagd14, mins, fmagd14, sc);


fmagu15: fmag;
fmagd15: fmag;
dmagu15: dmag;
dmagd15: dmag;
rf03: line=(drrf, rfc, drrf);
rf04: line=(drrf, rfc, drrf);

cell15: line = (sa, cpshort, sb, fmagu15, mins, dmagu15, drifta, rf03, driftb, rf04, dmidls, cplong, drifte, dmagd15, mins, fmagd15, sc);

fmagu16: fmag;
fmagd16: fmag;
dmagu16: dmag;
dmagd16: dmag;
rf05: line=(drrf, rfc, drrf);
rf06: line=(drrf, rfc, drrf);

cell16: line = (sa, cpshort, sb, fmagu16, mins, dmagu16, drifta, rf05, driftb, rf06, dmidls, cplong, drifte, dmagd16, mins, fmagd16, sc);

fmagu17: fmag;
fmagd17: fmag;
dmagu17: dmag;
dmagd17: dmag;
rf07: line=(drrf, rfc, drrf);
rf08: line=(drrf, rfc, drrf);

cell17: line = (sa, cpshort, sb, fmagu17, mins, dmagu17, drifta, rf07, driftb, rf08, dmidls, cplong, drifte, dmagd17, mins, fmagd17, sc);

fmagu18: fmag;
fmagd18: fmag;
dmagu18: dmag;
dmagd18: dmag;
rf09: line=(drrf, rfc, drrf);
rf10: line=(drrf, rfc, drrf);

cell18: line = (sa, cpshort, sb, fmagu18, mins, dmagu18, drifta, rf09, driftb, rf10, dmidls, cplong, drifte, dmagd18, mins, fmagd18, sc);

fmagu19: fmag;
fmagd19: fmag;
dmagu19: dmag;
dmagd19: dmag;
rf11: line=(drrf, rfc, drrf);
rf12: line=(drrf, rfc, drrf);

cell19: line = (sa, cpshort, sb, fmagu19, mins, dmagu19, drifta, rf11, driftb, rf12, dmidls, cplong, drifte, dmagd19, mins, fmagd19, sc);

fmagu20: fmag;
fmagd20: fmag;
dmagu20: dmag;
dmagd20: dmag;
rf13: line=(drrf, rfc, drrf);
rf14: line=(drrf, rfc, drrf);

cell20: line = (sa, cpshort, sb, fmagu20, mins, dmagu20, drifta, rf13, driftb, rf14, dmidls, cplong, drifte, dmagd20, mins, fmagd20, sc);

fmagu21: fmag;
fmagd21: fmag;
dmagu21: dmag;
dmagd21: dmag;
rf15: line=(drrf, rfc, drrf);
rf16: line=(drrf, rfc, drrf);

cell21: line = (sa, cpshort, sb, fmagu21, mins, dmagu21, drifta, rf15, driftb, rf16, dmidls, cplong, drifte, dmagd21, mins, fmagd21, sc);

fmagu22: fmag;
fmagd22: fmag;
dmagu22: dmag;
dmagd22: dmag;
rf17: line=(drrf, rfc, drrf);
rf18: line=(drrf, rfc, drrf);

cell22: line = (sa, cpshort, sb, fmagu22, mins, dmagu22, drifta, rf17, driftb, rf18, dmidls, cplong, drifte, dmagd22, mins, fmagd22, sc);

fmagu23: fmag;
fmagd23: fmag;
dmagu23: dmag;
dmagd23: dmag;
rf19: line=(drrf, rfc, drrf);
rf20: line=(drrf, rfc, drrf);

cell23: line = (sa, cpshort, sb, fmagu23, mins, dmagu23, drifta, rf19, driftb, rf20, dmidls, cplong, drifte, dmagd23, mins, fmagd23, sc);

fmagu24: fmag;
fmagd24: fmag;
dmagu24: dmag;
dmagd24: dmag;
rf21: line=(drrf, rfc, drrf);
rf22: line=(drrf, rfc, drrf);

cell24: line = (sa, cpshort, sb, fmagu24, mins, dmagu24, drifta, rf21, driftb, rf22, dmidls, cplong, drifte, dmagd24, mins, fmagd24, sc);

booster: line=(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24);

beam, particle=proton, energy=0.8+pmass;
)foo");

    MadX_reader reader;
    reader.parse(booster_madx);
    Lattice lattice(reader.get_lattice("booster"));

    Lattice_simulator::tune_circular_lattice(lattice);
    return lattice;
}

Propagator
create_propagator(Lattice lattice)
{
    Independent_stepper_elements stepper(1);
    Propagator prop(lattice, stepper);
    return prop;
}

TEST_CASE("create_propagator")
{
    Lattice lattice(get_lattice());
    Propagator p(create_propagator(lattice));
}

// create and initialize the bunch simulator for propagation
Bunch_simulator
create_simulator(Lattice const& lattice)
{
  
  // create the simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(
							      lattice.get_reference_particle(),
							      macroparticles,
							      realparticles,
							      Commxx(),
							      spectparticles);

    // initializate particles in the simulator
    auto& bunch = sim.get_bunch(0, 0);

    auto bp = bunch.get_local_particles(ParticleGroup::regular);
    for (int i=0; i<macroparticles; ++i) {
      for (int j=0; j<6; ++j) {
	bp(i, j) = 0.0;
      }
    }
    auto sbp = bunch.get_local_particles(ParticleGroup::spectator);
    for (int i=0; i<spectparticles; ++i) {
      for (int j=0; j<6; ++j) {
	sbp(i, j) = 0.0;
      }
    }

    bunch.checkin_particles();
    
    return sim;

}


void
test_particles(Bunch_simulator& sim, Lattice const& lattice, int turn, int step)
{
  // test that particle 0 from regular and spectator particles matches
  auto& bunch = sim.get_bunch(0, 0);
  auto bp = bunch.get_local_particles(ParticleGroup::regular);
  auto sp = bunch.get_local_particles(ParticleGroup::spectator);

  auto const& lems = lattice.get_elements();
  // steps start at 1 for some reason

  auto it = lems.begin();
  std::advance(it, step-1);

#if 0
  std::cout << "step: " << step << std::endl;
  std::cout << it->get_name() << ": " << it->get_type_name() << std::endl;
  std::cout << "0: " << std::scientific << std::setprecision(16) << bp(0, 0) << " <-> " << sp(0, 0) << std::endl;
  std::cout << "1: " << std::scientific << std::setprecision(16) << bp(0, 1) << " <-> " << sp(0, 1) << std::endl;
  std::cout << "2: " << std::scientific << std::setprecision(16) << bp(0, 2) << " <-> " << sp(0, 2) << std::endl;
  std::cout << "3: " << std::scientific << std::setprecision(16) << bp(0, 3) << " <-> " << sp(0, 3) << std::endl;
  std::cout << "4: " << std::scientific << std::setprecision(16) << bp(0, 4) << " <-> " << sp(0, 4) << std::endl;
  std::cout << "5: " << std::scientific << std::setprecision(16) << bp(0, 5) << " <-> " << sp(0, 5) << std::endl;
  std::cout << std::endl;
#endif

  REQUIRE(bp(0, 0) == sp(0, 0));
  REQUIRE(bp(0, 1) == sp(0, 1));
  REQUIRE(bp(0, 2) == sp(0, 2));
  REQUIRE(bp(0, 3) == sp(0, 3));
  REQUIRE(bp(0, 4) == sp(0, 4));
  REQUIRE(bp(0, 5) == sp(0, 5));
}

TEST_CASE("propagate_particles")
{
  Lattice lattice(get_lattice());

    auto sim = create_simulator(lattice);

    sim.reg_prop_action_step_end(test_particles);

    Propagator p(create_propagator(lattice));

    Logger simlog(0, LoggerV::INFO_TURN);

    p.propagate(sim, simlog, nturns);
}

