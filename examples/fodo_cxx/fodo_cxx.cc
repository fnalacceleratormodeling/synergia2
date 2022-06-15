// standalone c++ version of the fodo.py example

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/populate.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/simulation/checkpoint.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include <iomanip>
#include <iostream>
#include <string>

#include "synergia/collective/space_charge_3d_open_hockney.h"

#include "fodo_cxx_options.h"

Lattice
get_lattice()
{
  static std::string fodo_madx(R"foo(
beam, particle=proton,pc=3.0;

o: drift, l=8.0;
f: quadrupole, l=2.0, k1=0.071428571428571425;
d: quadrupole, l=2.0, k1=-0.071428571428571425;

fodo: sequence, l=20.0, refer=entry;
fodo_1: f, at=0.0;
fodo_2: o, at=2.0;
fodo_3: d, at=10.0;
fodo_4: o, at=12.0;
endsequence;
)foo");

  MadX_reader reader;
  reader.parse(fodo_madx);
  return reader.get_lattice("fodo");
}

void
print_bunch_statistics(Bunch const& bunch, Logger& logger)
{
  karray1d bunch_means(Core_diagnostics::calculate_mean(bunch));
  karray1d_row bunch_stds(Core_diagnostics::calculate_std(bunch, bunch_means));

  logger << "bunch means" << std::endl;
  for (int i = 0; i < 6; ++i) {
    logger << i << ":  " << std::setprecision(15) << bunch_means(i)
           << std::endl;
  }
  logger << "bunch stds" << std::endl;
  for (int i = 0; i < 6; ++i) {
    logger << i << ":  " << std::setprecision(15) << bunch_stds(i) << std::endl;
  }
}

int
run(Fodo_cxx_options opts)
{
  int gridx = opts.gridx;
  int gridy = opts.gridy;
  int gridz = opts.gridz;
  double real_particles = opts.real_particles;
  int macroparticles = opts.macroparticles;
  int turns = opts.turns;

  Logger screen(0, LoggerV::INFO);
  // Logger screen(0, LoggerV::DEBUG);

  screen << "gridx: " << gridx << std::endl;
  screen << "gridy: " << gridy << std::endl;
  screen << "gridz: " << gridz << std::endl;
  screen << "macroparticles: " << macroparticles << std::endl;
  screen << "real_particles: " << real_particles << std::endl;

  Lattice lattice = get_lattice();
  screen << "Read lattice, length: " << lattice.get_length() << ", "
         << lattice.get_elements().size() << " elements" << std::endl;

  Four_momentum four_momentum(pconstants::mp);
  four_momentum.set_momentum(3.0);
  Reference_particle refpart(1, four_momentum);
  lattice.set_reference_particle(refpart);

  // space charge
  Space_charge_3d_open_hockney_options sc_ops(gridx, gridy, gridz);
  sc_ops.comm_group_size = 1;

  // stepper
  Split_operator_stepper_elements stepper(sc_ops, 1);

  // Propagator
  Propagator propagator(lattice, stepper);
  // propagator.print_steps(screen);

  // print slices
  for (auto const& slice : propagator.get_lattice_element_slices())
    screen << slice.as_string() << "\n";

  // bunch simulator
  auto sim = Bunch_simulator::create_single_bunch_simulator(
    lattice.get_reference_particle(), macroparticles, real_particles, Commxx());

  karray1d means("means", 6);
  for (int i = 0; i < 6; ++i) means(i) = 0.0;

  static const double covariance_matrix[6][6] = {
    {3.0509743977035345e-05, 2.2014134466660509e-06, 0, 0, 0, 0},
    {2.2014134466660509e-06, 1.9161816525115869e-07, 0, 0, 0, 0},
    {0, 0, 7.5506914064526925e-06, -6.6846812465678249e-07, 0, 0},
    {0, 0, -6.6846812465678249e-07, 1.9161816525115867e-07, 0, 0},
    {0, 0, 0, 0, 0.00016427607645871527, 0},
    {0, 0, 0, 0, 0, 1e-08}};

  karray2d_row covariances("covariances", 6, 6);
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) covariances(i, j) = covariance_matrix[i][j];

  Random_distribution dist(1234567, Commxx());

  auto& bunch = sim.get_bunch();
  populate_6d(dist, bunch, means, covariances);

  screen << "Statistics before propagation" << std::endl;

  bunch.checkout_particles();
  print_bunch_statistics(bunch, screen);

  // diagnostics
  /*
  Diagnostics_bulk_track diag_bulk_track(
    "tracks.h5", (100 < macroparticles) ? 100 : macroparticles, 0);
  sim.reg_diag_per_turn(diag_bulk_track);
  */

  Diagnostics_particles diag_particles("particles.h5", macroparticles);
  sim.reg_diag_per_turn(diag_particles);

  /*
  Diagnostics_full2 diag_full2("diag.h5");
  sim.reg_diag_per_turn(diag_full2);

  screen << "Statistics before propagation" << std::endl;
  */

  bunch.checkout_particles();
  print_bunch_statistics(bunch, screen);
  // propagate
  Logger proplogger = Logger(0, LoggerV::INFO_TURN);
  propagator.propagate(sim, proplogger, turns);

  bunch.checkout_particles();
  screen << "Statistics after propagate" << std::endl;
  print_bunch_statistics(bunch, screen);

  syn::checkpoint_save(propagator, sim);

  return 0;
}

int
main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);

  Fodo_cxx_options opts(argc, argv);

  // layout_test();
  run(opts);

#ifdef SIMPLE_TIMER
  Logger logger(0);
  simple_timer_print(logger);
#endif

  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
