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


constexpr double toffs = 1.0e-3;
constexpr double dpoffs = 1.0e-4;
constexpr double cdtoffs =  0.05;
constexpr double dpopoffs = 2.2e-4;

Lattice
get_lattice(std::string const& lattice_madx)
{

    MadX_reader reader;
    reader.parse(lattice_madx);
    return reader.get_lattice("element");
}

Propagator
create_propagator(Lattice lattice, int nsteps)
{
    Independent_stepper_elements stepper(nsteps);
    Propagator prop(lattice, stepper);
    return prop;
}

#if 0
TEST_CASE("create_propagator")
{
    Lattice lattice(get_lattice(iota_madx1));
    Propagator p1(create_propagator(lattice, 1));
    Propagator p2(create_propagator(lattice, 3));
}
#endif

void
fill_bunch(Bunch& bunch)
{
    auto bp = bunch.get_local_particles(ParticleGroup::regular);
    for (int i=0; i<macroparticles; ++i) {
      for (int j=0; j<6; ++j) {
	bp(i, j) = 0.0;
      }
    }
    // fill in test particles
    // particle 0 remains at 0
    bp(1, 0) = toffs; // particle 1 offset in x
    bp(2, 1) = dpoffs; // particle 2 momentum offset
    bp(3, 2) = toffs;
    bp(4, 3) = dpoffs;
    bp(5, 4) = cdtoffs;
    bp(6, 5) = dpopoffs;

    bp(7, 0) = -toffs;
    bp(7, 5) = dpopoffs;

    bp(8, 1) = -dpoffs;
    bp(8, 5) = dpopoffs;

    bp(9, 2) = toffs;
    bp(9, 3) = -dpoffs;
    bp(9, 5) = -dpopoffs;

    bp(10, 0) = -toffs;
    bp(10, 1) = -dpoffs;
    bp(10, 2) = toffs;
    bp(10, 3) = dpoffs;
    bp(10, 5) = dpopoffs;

    bp(11, 0) = toffs;
    bp(11, 1) = dpoffs;
    bp(11, 2) = -toffs;
    bp(11, 3) = -dpoffs;
    bp(11, 5) = -dpopoffs;
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
							      Commxx());

    // initializate particles in the simulator
    auto& bunch = sim.get_bunch(0, 0);

    fill_bunch(bunch);
    bunch.checkin_particles();
    
    return sim;

}


#if 0
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
#endif

void
propagate_particles(Lattice lattice, Bunch_simulator& sim, int nsteps)
{
// sim.reg_prop_action_step_end(test_particles);

  Propagator p(create_propagator(lattice, nsteps));

  Logger simlog(0, LoggerV::INFO_STEP);

  p.propagate(sim, simlog, 1);
}

void
run_test(std::string lattice_madx)
{
  Lattice lattice(get_lattice(lattice_madx));
  auto sim1 = create_simulator(lattice);
  auto sim2 = create_simulator(lattice);
  auto sim3 = create_simulator(lattice);
  auto sim5 = create_simulator(lattice);

  Propagator p1(create_propagator(lattice, 1));
  Propagator p2(create_propagator(lattice, 2));
  Propagator p3(create_propagator(lattice, 3));
  Propagator p5(create_propagator(lattice, 5));

  Logger simlog(0, LoggerV::INFO_STEP);
  p1.propagate(sim1, simlog, 1);
  p2.propagate(sim2, simlog, 1);
  p3.propagate(sim3, simlog, 1);
  p5.propagate(sim5, simlog, 1);
  
  auto& bunch1 = sim1.get_bunch(0, 0);
  bunch1.checkout_particles();
  auto& bunch2 = sim2.get_bunch(0, 0);
  bunch2.checkout_particles();
  auto& bunch3 = sim3.get_bunch(0, 0);
  bunch3.checkout_particles();
  auto& bunch5 = sim5.get_bunch(0, 0);
  bunch5.checkout_particles();

  auto bp1 = bunch1.get_local_particles();
  auto bp2 = bunch2.get_local_particles();
  auto bp3 = bunch3.get_local_particles();
  auto bp5 = bunch5.get_local_particles();

  for (int ip=0; ip<11; ++ip) {
    for (int j=0; j<6; ++j) {

      std::cout << "steps 1-2  (" << ip << "," << j << "): " << std::scientific << std::setprecision(16) << bp1(ip, j) << " <-> " << bp2(ip, j) << std::endl;
      if((abs(bp1(ip, j)) < 1.0e-9) && (abs(bp2(ip, j) < 1.0e-9))) {
	REQUIRE_THAT(bp1(ip, j), Catch::Matchers::WithinAbs(bp2(ip, j), 1.0e-14));
      } else {
	REQUIRE_THAT(bp1(ip, j), Catch::Matchers::WithinRel(bp2(ip, j), 1.0e-14));
      }

      std::cout << "steps 1-3  (" << ip << "," << j << "): " << std::scientific << std::setprecision(16) << bp1(ip, j) << " <-> " << bp3(ip, j) << std::endl;
      if((abs(bp1(ip, j)) < 1.0e-9) && (abs(bp3(ip, j) < 1.0e-9))) {
	REQUIRE_THAT(bp1(ip, j), Catch::Matchers::WithinAbs( bp3(ip, j), 1.0e-14));
      } else {
	REQUIRE_THAT(bp1(ip, j), Catch::Matchers::WithinRel( bp3(ip, j), 1.0e-14));
      }

      std::cout << "steps 1-5  (" << ip << "," << j << "): " << std::scientific << std::setprecision(16) << bp1(ip, j) << " <-> " << bp5(ip, j) << std::endl;
      if((abs(bp1(ip, j)) < 1.0e-9) && (abs(bp5(ip, j) < 1.0e-9))) {
	REQUIRE_THAT(bp1(ip, j), Catch::Matchers::WithinAbs(bp5(ip, j), 1.0e-14));
      } else {
	REQUIRE_THAT(bp1(ip, j), Catch::Matchers::WithinRel(bp5(ip, j), 1.0e-14));
      }
    }
  }
}

// this one better work
TEST_CASE("drift")
{
  static std::string lattice_madx(R"foo(
d: drift, l=2.0;
element: sequence, l=2.0, refer=entry;
 d, at=0.0;
 endsequence;
beam, particle=proton, energy = 0.00250+pmass;
)foo");

  run_test(lattice_madx);

}

TEST_CASE("iota_sbend")
{
  static std::string lattice_madx(R"foo(
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
element: sequence, l=0.3911403725, refer=entry;
m1r, at=0.0;
 endsequence;
beam, particle=proton, energy = 0.00250+pmass;
)foo");

  run_test(lattice_madx);

}


