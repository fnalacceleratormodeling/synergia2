#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_particles.h"

#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/bunch_simulator.h"

#include "synergia/utils/utils.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"

int macroparticles = 16;
int spectparticles = 4;
double realparticles = 285452129459.3449; // 0.1 mA in IOTA
int nturns = 5;

Lattice
get_lattice()
{
    static std::string iota_madx(R"foo(
kqa1r := kq01;
kq01 = -7.72652301;
kqa2r := kq02;
kq02 = 12.28401222;
kqa3r := kq03;
kq03 = -12.43016989;
kqa4r := kq04;
kq04 = 20.16074347;
kqb1r := kq05;
kq05 = -10.24365752;
ibpm: monitor;
ibpma3r: ibpm;
mseqare: marker;
mphm1ri: marker;
sqa2r: quadrupole,l:= 0.1,k1s:= 0;
qa4r: quadrupole,l:= 0.21,k1:=kqa4r ;
dedge30: sextupole, l=0.0, k2=0.0;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;



!iota: sequence, l = 39.95567226;
iota: sequence, l = 3.508540575;
!mseqari, at = 0;
!ibpma1, at = 0.02;
!qa1r, at = 1.0175;
!qa2r, at = 1.3625;
!sqa1r, at = 1.99;
!ibpma2r, at = 2.095;
!qa3r, at = 2.2975;
qa4r, at = 2.6525; ! this one has to be there for the error to occur
ibpma3r, at = 2.865;
sqa2r, at = 2.97;
mseqare, at = 3.0405;
mphm1ri, at = 3.0405;
!dedge30, at = 3.117400202;
m1r, at = 3.312970389;
endsequence;
beam, particle=proton, energy = 0.00250 + pmass;
)foo");

    MadX_reader reader;
    reader.parse(iota_madx);
    return reader.get_lattice("iota");
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

TEST_CASE("propagate_particles")
{
  Lattice lattice(get_lattice());

    auto sim = create_simulator(lattice);

    sim.reg_prop_action_step_end(test_particles);

    Propagator p(create_propagator(lattice));

    Logger simlog(0, LoggerV::INFO_STEP);

    p.propagate(sim, simlog, nturns);
}

