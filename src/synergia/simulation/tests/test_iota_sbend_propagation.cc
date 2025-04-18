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
int spectparticles = 16;
double realparticles = 285452129459.3449; // 0.1 mA in IOTA
int nturns = 5;

static std::string iota_madx1(R"foo(
mphm1ri: marker;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;

iota: sequence, l = 0.3911403725;
mphm1ri, at = 0.0;
m1r, at = 0.3911403725/2;
endsequence;
beam, particle=proton, energy = 0.00250 + pmass;
)foo");

static std::string iota_madx2(R"foo(
mphm1ri: marker;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756, k1=0.0;

iota: sequence, l = 0.3911403725;
mphm1ri, at = 0.0;
m1r, at = 0.3911403725/2;
endsequence;
beam, particle=proton, energy = 0.00250 + pmass;
)foo");

static std::string iota_madx3(R"foo(
mphm1ri: marker;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756, e1=0.5235987756/2, e2=0.5235987756/2;

iota: sequence, l = 0.3911403725;
mphm1ri, at = 0.0;
m1r, at = 0.3911403725/2;
endsequence;
beam, particle=proton, energy = 0.00250 + pmass;
)foo");

Lattice
get_lattice(std::string const& lattice_madx)
{

    MadX_reader reader;
    reader.parse(lattice_madx);
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
    Lattice lattice(get_lattice(iota_madx1));
    Propagator p(create_propagator(lattice));
}

// create and initialize the bunch simulator for propagation
Bunch_simulator
create_simulator(Lattice const& lattice, double cdt0 = 0.0)
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
      // bp(0, 4) = -8.8817841970012523e-16;
       bp(0, 4) = cdt0;
    }
    auto sbp = bunch.get_local_particles(ParticleGroup::spectator);
    for (int i=0; i<spectparticles; ++i) {
      for (int j=0; j<6; ++j) {
	sbp(i, j) = 0.0;
      }
      // sbp(0, 4) = -8.8817841970012523e-16;
      sbp(0, 4) = cdt0;
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

TEST_CASE("propagate_particles")
{
  Lattice lattice(get_lattice(iota_madx1));

    auto sim = create_simulator(lattice);

    sim.reg_prop_action_step_end(test_particles);

    Propagator p(create_propagator(lattice));

    Logger simlog(0, LoggerV::INFO_STEP);

    p.propagate(sim, simlog, nturns);
}

TEST_CASE("propagate_particles1")
{
  Lattice lattice(get_lattice(iota_madx2));

  auto sim = create_simulator(lattice, -8.8817841970012523e-16);

  sim.reg_prop_action_step_end(test_particles);

  Propagator p(create_propagator(lattice));

  Logger simlog(0, LoggerV::INFO_STEP);

  p.propagate(sim, simlog, nturns);
}

TEST_CASE("propagate_particles2")
{
  Lattice lattice(get_lattice(iota_madx3));

  auto sim = create_simulator(lattice, -8.8817841970012523e-16);

  sim.reg_prop_action_step_end(test_particles);

  Propagator p(create_propagator(lattice));

  Logger simlog(0, LoggerV::INFO_STEP);

  p.propagate(sim, simlog, nturns);
}
