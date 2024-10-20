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
! File fobodobo_s.lat
! 
! Written for the January, 2007 USPAS.
! To be used in conjunction with CHEF.
! 
! Send complaints to the author: Leo Michelotti
! 
! Add an RF cavity (EGS) 01/30/2009
!
! ------------------
! Parameters
! ------------------
n           :=   32;                   !         : number of cells
bendangle   := twopi/(2*n);           ! [rad]   : dipole bend angle
focus       :=   7;                   ! [m]     : focal length of equivalent 
                                     !         :   thin quad
sepn        :=  10;                   ! [m]     : distance between quad centers
quadlength  :=   0.2;                 ! [m]     : quadrupole length
strength    := 1/(focus*quadlength);  ! [m**-2] : quadrupole strength
                                     !         :   = B'/brho, where
                                     !         :   brho = momentum/0.299792458
pct         :=   0.4;                 !         : fraction of space between
                                     !         :   quads occupied by dipole
bendlength  := pct*(sepn-quadlength); ! [m]     : length of dipole
! bendlength := 0.5*(10-2.0) = 4.0
driftlength := (sepn-quadlength-bendlength)/2;
! driftlenth := (10-2.0-4.0)/2 = 2.0
! harmonic number = 80  harmonic number, yields 2 meter wavelength
! the actual frequence is harmno following
harmno:=32;
lambda = (n*2*sepn)/harmno;

!hvoltage = 12500
hvoltage = 0.000001;

cavlen = 0.0; ! rf cavity length 1 meter, (half bucket length)
shortdlen = (driftlength - cavlen)/2.0;   ! 0.97 m
! this lattice has 32 copies of a 20 m cell.

! ------------------
! Elements
! ------------------

o: drift,      l=driftlength;
os: drift,      l=shortdlen;
f: quadrupole, l=quadlength, k1=strength;
d: quadrupole, l=quadlength, k1=(-strength);
b: sbend,      l=bendlength, angle=bendangle; ! this one provokes the error with spectators
!b: sbend,      l=bendlength, angle=bendangle, e1=bendangle/2, e2=bendangle/2; ! this one doesn't get an error
r: rfcavity,l=cavlen, volt=hvoltage, harmon=harmno, lag=0.0;

! chromaticity adjusters
sf: sextupole, l=0, k2=0.0;
sd: sextupole, l=0, k2=-0.0;


! ------------------
! Lattices
! ------------------
fobodobo:  line=( f, sf, o, b, o, d, sd, o, b, o );
fobrdobo:  line=( f, sf,  o, b, os, r, os, d, sd, o, b, o);
model:     line=( fobrdobo,31*fobodobo );

beam, particle=proton, energy=0.8+pmass;
)foo");

    MadX_reader reader;
    reader.parse(booster_madx);
    Lattice lattice(reader.get_lattice("model"));

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
create_simulator(Lattice const& lattice, int spart=0)
{
  
  // create the simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(
							      lattice.get_reference_particle(),
							      macroparticles,
							      realparticles,
							      Commxx(),
							      spart);

    // initializate particles in the simulator
    auto& bunch = sim.get_bunch(0, 0);

    auto bp = bunch.get_local_particles(ParticleGroup::regular);
    for (int i=0; i<macroparticles; ++i) {
      for (int j=0; j<6; ++j) {
	bp(i, j) = 0.0;
      }
    }
    auto sbp = bunch.get_local_particles(ParticleGroup::spectator);
    for (int i=0; i<bunch.get_local_num(ParticleGroup::spectator); ++i) {
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

  std::cout << "step: " << step << std::endl;
  std::cout << it->get_name() << ": " << it->get_type_name() << std::endl;
  std::cout << "0: " << std::scientific << std::setprecision(16) << bp(0, 0) << " <-> " << sp(0, 0) << std::endl;
  std::cout << "1: " << std::scientific << std::setprecision(16) << bp(0, 1) << " <-> " << sp(0, 1) << std::endl;
  std::cout << "2: " << std::scientific << std::setprecision(16) << bp(0, 2) << " <-> " << sp(0, 2) << std::endl;
  std::cout << "3: " << std::scientific << std::setprecision(16) << bp(0, 3) << " <-> " << sp(0, 3) << std::endl;
  std::cout << "4: " << std::scientific << std::setprecision(16) << bp(0, 4) << " <-> " << sp(0, 4) << std::endl;
  std::cout << "5: " << std::scientific << std::setprecision(16) << bp(0, 5) << " <-> " << sp(0, 5) << std::endl;
  std::cout << std::endl;

  REQUIRE(bp(0, 0) == sp(0, 0));
  REQUIRE(bp(0, 1) == sp(0, 1));
  REQUIRE(bp(0, 2) == sp(0, 2));
  REQUIRE(bp(0, 3) == sp(0, 3));
  REQUIRE(bp(0, 4) == sp(0, 4));
  REQUIRE(bp(0, 5) == sp(0, 5));
}

TEST_CASE("propagate_particles_with_spectators")
{
  Lattice lattice(get_lattice());

    auto sim = create_simulator(lattice, spectparticles);

    sim.reg_prop_action_step_end(test_particles);

    Propagator p(create_propagator(lattice));

    Logger simlog(0, LoggerV::INFO_TURN);

    p.propagate(sim, simlog, nturns);
}

std::list<double> bunch_cdts;

void
save_cdt(Bunch_simulator& sim, Lattice const& lattice, int turn, int step)
{
  // test that particle 0 from regular and spectator particles matches
  auto& bunch = sim.get_bunch(0, 0);
  auto bp = bunch.get_local_particles(ParticleGroup::regular);
  auto sp = bunch.get_local_particles(ParticleGroup::spectator);

  bunch_cdts.push_back(bp(0, 4));
}

 TEST_CASE("propagate_particles_without_spectators")
{
  Lattice lattice(get_lattice());

  auto sim = create_simulator(lattice, 0);
  sim.reg_prop_action_step_end(save_cdt);

  Propagator p(create_propagator(lattice));

  Logger simlog(0, LoggerV::INFO_TURN);

  p.propagate(sim, simlog, nturns);

  std::list<double> cdt_nospect(bunch_cdts);
  bunch_cdts.clear();
  auto sim2 = create_simulator(lattice, spectparticles);
  sim2.reg_prop_action_step_end(save_cdt);

  p.propagate(sim2, simlog, nturns);

  std::list<double> cdt_yesspect(bunch_cdts);

  std::cout << "cdt_nospect.size(): " << cdt_nospect.size() << std::endl;
  std::cout << "cdt_yesspect.size(): " << cdt_yesspect.size() << std::endl;

  REQUIRE(cdt_nospect.size() == cdt_yesspect.size());

  std::list<double>::const_iterator it1=cdt_nospect.begin();
  std::list<double>::const_iterator it2=cdt_yesspect.begin();
  int stpcnt = 1;
  for(	; it1!=cdt_nospect.end(); ++it1, ++it2, ++stpcnt) {
    if (*it1 != *it2) {
      std::cout << "step: " << stpcnt << " mismatch: " <<
	  std::scientific << std::setprecision(16) << *it1 << " <-> " << *it2 << std::endl;
    }
    REQUIRE(*it1 == *it2);
  }
  
}
   
TEST_CASE("long_term_stability")
{
  Lattice lattice(get_lattice());

  auto sim = create_simulator(lattice, 0);

  Propagator p(create_propagator(lattice));

  Logger simlog(0, LoggerV::INFO_TURN);

  p.propagate(sim, simlog, 1000);

  // test that particle 0 from regular and spectator particles matches
  auto& bunch = sim.get_bunch(0, 0);
  auto bp = bunch.get_local_particles(ParticleGroup::regular);

    std::cout << "0: " << std::scientific << std::setprecision(16) << bp(0, 0)  << std::endl;
    std::cout << "1: " << std::scientific << std::setprecision(16) << bp(0, 1) << " <-> " << std::endl;
  std::cout << "2: " << std::scientific << std::setprecision(16) << bp(0, 2) << " <-> " << std::endl;
  std::cout << "3: " << std::scientific << std::setprecision(16) << bp(0, 3) << " <-> " <<  std::endl;
  std::cout << "4: " << std::scientific << std::setprecision(16) << bp(0, 4) << " <-> " <<  std::endl;
  std::cout << "5: " << std::scientific << std::setprecision(16) << bp(0, 5) << " <-> " <<  std::endl;
  std::cout << std::endl;

}

TEST_CASE("long_term_stability1")
{
  Lattice lattice(get_lattice());

  auto sim = create_simulator(lattice, 4);

  Propagator p(create_propagator(lattice));

  Logger simlog(0, LoggerV::INFO_TURN);

  p.propagate(sim, simlog, 1000);

  // test that particle 0 from regular and spectator particles matches
  auto& bunch = sim.get_bunch(0, 0);
  auto bp = bunch.get_local_particles(ParticleGroup::regular);

    std::cout << "0: " << std::scientific << std::setprecision(16) << bp(0, 0)  << std::endl;
    std::cout << "1: " << std::scientific << std::setprecision(16) << bp(0, 1) << " <-> " << std::endl;
  std::cout << "2: " << std::scientific << std::setprecision(16) << bp(0, 2) << " <-> " << std::endl;
  std::cout << "3: " << std::scientific << std::setprecision(16) << bp(0, 3) << " <-> " <<  std::endl;
  std::cout << "4: " << std::scientific << std::setprecision(16) << bp(0, 4) << " <-> " <<  std::endl;
  std::cout << "5: " << std::scientific << std::setprecision(16) << bp(0, 5) << " <-> " <<  std::endl;
  std::cout << std::endl;

}
