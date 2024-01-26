#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/checkpoint.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/lsexpr.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/utils.h"
#include <array>
#include <iomanip>
#include <iostream>
#include <string>

#include "booster_fd_options.h"

#ifdef BUILD_FD_SPACE_CHARGE_SOLVER
#include "synergia/collective/space_charge_3d_fd.h"
#else
#include "synergia/collective/space_charge_rectangular.h"
#endif

Lattice
get_lattice()
{
    const std::string lattice_filename = "orbump_rampdown_0099.lsx";
    return Lattice(read_lsexpr_file(lattice_filename));
}

void
print_bunch_statistics(Bunch const& bunch, Logger& logger)
{
    karray1d bunch_means(Core_diagnostics::calculate_mean(bunch));
    karray1d_row bunch_stds(
        Core_diagnostics::calculate_std(bunch, bunch_means));

    logger << "bunch means" << std::endl;
    for (int i = 0; i < 6; ++i) {
        logger << i << ":  " << std::setprecision(15) << bunch_means(i)
               << std::endl;
    }
    logger << "bunch stds" << std::endl;
    for (int i = 0; i < 6; ++i) {
        logger << i << ":  " << std::setprecision(15) << bunch_stds(i)
               << std::endl;
    }
}

int
run(Booster_fd_options opts)
{
    int gridx = opts.gridx;
    int gridy = opts.gridy;
    int gridz = opts.gridz;
    double real_particles = opts.real_particles;
    int num_particles = opts.num_particles;
    int turns = opts.turns;
    double pipesizex = opts.pipesizex;
    double pipesizey = opts.pipesizey;
    double pipesizez = opts.pipesizez;

    Logger screen(0, LoggerV::INFO);
    // Logger screen(0, LoggerV::DEBUG);

    screen << "gridx: " << gridx << std::endl;
    screen << "gridy: " << gridy << std::endl;
    screen << "gridz: " << gridz << std::endl;
    screen << "real_particles: " << real_particles << std::endl;

    Lattice lattice = get_lattice();
    screen << "Read lattice, length: " << lattice.get_length() << ", "
           << lattice.get_elements().size() << " elements" << std::endl;

    auto refpart = lattice.get_reference_particle();
    auto beta = refpart.get_beta();

    // space charge
#ifdef BUILD_FD_SPACE_CHARGE_SOLVER
    Space_charge_3d_fd_options sc_ops(gridx, gridy, gridz);
    sc_ops.set_fixed_domain(
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{pipesizex, pipesizey, pipesizez / beta});
#else
    Space_charge_rectangular_options sc_ops(
        std::array<int, 3>{gridx, gridy, gridz},
        std::array<double, 3>{pipesizex, pipesizey, pipesizez});
#endif
    sc_ops.comm_group_size = 4;

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
        lattice.get_reference_particle(),
        num_particles,
        real_particles,
        Commxx());

    karray1d means("means", 6);
    for (int i = 0; i < 6; ++i)
        means(i) = 0.0;

    auto& bunch = sim.get_bunch();
    const std::string filename = "particles_booster_initial_onepercent.h5";
    bunch.read_file_legacy(filename);

    screen << "Statistics before propagation" << std::endl;

    bunch.checkout_particles();
    print_bunch_statistics(bunch, screen);

    // diagnostics

    Diagnostics_particles diag_particles("particles.h5", 1000);
    sim.reg_diag_per_turn(diag_particles);

    Diagnostics_full2 diag_full2("diag.h5");
    sim.reg_diag_per_turn(diag_full2);

    screen << "Statistics before propagation" << std::endl;

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
    synergia::initialize(argc, argv);

    Booster_fd_options opts;

    // layout_test();
    run(opts);

#ifdef SIMPLE_TIMER
    Logger logger(0);
    simple_timer_print(logger);
#endif

    synergia::finalize();
    return 0;
}
