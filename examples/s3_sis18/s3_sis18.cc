

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics_track.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/collective/space_charge_2d_open_hockney.h"

#include "synergia/utils/simple_timer.h"

using LV = LoggerV;

void print_statistics(Bunch & bunch, Logger & logger)
{
    // print particles after propagate
    bunch.checkout_particles();
    auto hparts = bunch.get_host_particles();

    double sum = 0;
    double mean[6] = {0, 0, 0, 0, 0, 0};
    double std[6] = {0, 0, 0, 0, 0, 0};

    for(int p=0; p<bunch.get_local_num(); ++p)
    {
        for(int i=0; i<6; ++i) 
        {
            sum += hparts(p, i);
            mean[i] += hparts(p, i);
        }
    }

    for (int i=0; i<6; ++i) mean[i]/=bunch.get_local_num();

    for(int p=0; p<bunch.get_local_num(); ++p)
    {
        for(int i=0; i<6; ++i) 
        {
            std[i] += (hparts(p, i) - mean[i]) * (hparts(p, i) - mean[i]);
        }
    }

    for (int i=0; i<6; ++i) std[i] = sqrt(std[i]/bunch.get_local_num());

    logger(LV::DEBUG) 
        << std::resetiosflags(std::ios::fixed)
        << std::setprecision(16)
        << "\n\nsum = " << sum << "\n"
        << std::setiosflags(std::ios::showpos | std::ios::scientific)
        ;


    for (int i=0; i<6; ++i) 
        logger(LV::DEBUG) << mean[i] << ", " << std[i] << "\n";

    logger(LV::DEBUG)
        << std::resetiosflags(std::ios::showpos | std::ios::scientific)
        << "\n";

    for (int p=0; p<4; ++p) bunch.print_particle(p, logger);

#if 0
    double g_sum = 0;
    MPI_Reduce(&sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, 0, bunch.get_comm());

    logger(LV::DEBUG) 
        << std::setprecision(8)
        << "\n\npropagated sum (reduced) = " << g_sum << "\n";
#endif

    logger << "\n";
}

int run()
{
    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_TURN);

    MadX_reader reader;
    auto lattice = reader.get_lattice("machine", "sis18.madx");
    Lattice_simulator::tune_circular_lattice(lattice);

    auto const & ref = lattice.get_reference_particle();

    //lattice.print(screen);

    screen(LV::INFO) 
        << "reference momentum = " << ref.get_momentum() << " GeV\n";

    // space charge
    Space_charge_2d_open_hockney_options sc_ops(64, 64, 64);
    sc_ops.comm_group_size = 1;

    // stepper
    //Independent_stepper_elements stepper(1);
    //Split_operator_stepper_elements stepper(1, sc_ops);
    Split_operator_stepper stepper(sc_ops, 71);

    // Propagator
    Propagator propagator(lattice, stepper);
    //propagator.print_steps(screen);

    // bunch simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            lattice.get_reference_particle(), 1024 * 1024 * 1, 2.94e10,
            Commxx() );

#if 0
    // propagate actions
    double cdt = 0.0;
    sim.reg_prop_action_step_end( [&cdt](Bunch_simulator& sim, Lattice&, int, int step, void*) { 
        cdt += sim.get_bunch().get_design_reference_particle().get_state()[4]; 
    }, nullptr );
#endif

    // propagate options
    sim.set_turns(0, 2); // (start, num_turns)

    auto & bunch = sim.get_bunch();

#if 0
    // populate particle data
    karray1d means("means", 6);
    for (int i=0; i<6; ++i) means(i) = 0.0;

    karray2d covariances("covariances", 6, 6);
    for (int i=0; i<6; ++i)
        for (int j=0; j<6; ++j)
            covariances(i, j) = 0.0;

    covariances(0,0) = 1e-2;
    covariances(1,1) = 1e-2;
    covariances(2,2) = 1e-2;
    covariances(3,3) = 1e-2;
    covariances(4,4) = 1e-2;
    covariances(5,5) = 1e-2;

    Random_distribution dist(5, Commxx());
    populate_6d(dist, bunch, means, covariances);
#endif

#if 1
    // or read from file
    bunch.read_file("turn_particles_0000.h5");
#endif

    // statistics before propagate
    print_statistics(bunch, screen);


#if 0
    // diagnostics
    Diagnostics_track diag_track(2, "part_2_track.h5");
    sim.reg_diag_per_turn("track_2", diag_track);

    Diagnostics_bulk_track diag_bulk_track(6, 0, "bulk_track.h5");
    sim.reg_diag_per_turn("bulk_track", diag_bulk_track);
#endif

    // propagate
    propagator.propagate(sim, simlog);

    // statistics after propagate
    print_statistics(bunch, screen);
    simple_timer_print(screen);

    return 0;
}


int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    run();

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

