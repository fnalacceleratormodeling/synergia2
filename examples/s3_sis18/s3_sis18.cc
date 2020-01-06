
#include "synergia/foundation/physical_constants.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/bunch/populate.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/bunch/diagnostics_loss.h"
//#include "synergia/bunch/diagnostics_track.h"
//#include "synergia/bunch/diagnostics_bulk_track.h"
//#include "synergia/bunch/diagnostics_full2.h"

#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/lsexpr.h"

#include "synergia/collective/space_charge_2d_open_hockney.h"

#include "synergia/utils/simple_timer.h"

using LV = LoggerV;

void print_statistics(Bunch & bunch, Logger & logger)
{
    // print particles after propagate
    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std  = Core_diagnostics::calculate_std(bunch, mean);

    logger(LV::DEBUG) 
        << std::resetiosflags(std::ios::fixed)
        << std::setprecision(16)
        << std::setiosflags(std::ios::showpos | std::ios::scientific)
        << "\n"
        ;

    for (int i=0; i<6; ++i) 
        logger(LV::DEBUG) << mean[i] << ", " << std[i] << "\n";

    logger(LV::DEBUG)
        << std::resetiosflags(std::ios::showpos | std::ios::scientific)
        << "\n";

    for (int p=0; p<4; ++p) bunch.print_particle(p, logger);

    logger << "\n";
}

int run()
{
    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_TURN);

#if 0
    MadX_reader reader;
    auto lattice = reader.get_lattice("machine", "sis18.madx");
#endif

    auto lsexpr = read_lsexpr_file("sis18-6.lsx");
    Lattice lattice(lsexpr);

    lattice.set_all_string_attribute("extractor_type", "libff");
    // lattice.print(screen);

    // tune the lattice
    Lattice_simulator::tune_circular_lattice(lattice);

    // get the reference particle
    auto const & ref = lattice.get_reference_particle();

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
            //lattice.get_reference_particle(), 1024 * 1024 * 4, 2.94e10,
            lattice.get_reference_particle(), 4194394, 2.94e10,
            Commxx() );

#if 0
    // propagate actions
    double cdt = 0.0;
    sim.reg_prop_action_step_end( [&cdt](Bunch_simulator& sim, Lattice&, int, int step, void*) { 
        cdt += sim.get_bunch().get_design_reference_particle().get_state()[4]; 
    }, nullptr );
#endif

    // propagate options
    sim.set_turns(0, 1); // (start, num_turns)

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
    bunch.read_file("turn_particles_0000_4M.h5");
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

#if 0
    Diagnostics_full2 diag_full2("diag_full.h5");
    sim.reg_diag_per_turn("full2", diag_full2);
    sim.reg_diag_per_turn("full3", Diagnostics_full2("diag_full3.h5"));
#endif

    // propagate
    propagator.propagate(sim, simlog);

    // statistics after propagate
    print_statistics(bunch, screen);
    simple_timer_print(screen);

    return 0;
}

std::string ar_name()
{
    std::stringstream ss;
    ss << "cp-" << Commxx().rank() << ".json";
    return ss.str();
}

void save()
{
    std::ofstream os(ar_name());
    cereal::JSONOutputArchive ar(os);

    auto sim = Bunch_simulator::create_single_bunch_simulator(
            Reference_particle(), 4194394, 2.94e10,
            Commxx() );

    Bunch& b = sim.get_bunch();

    Diagnostics_dummy diag;
    b.add_diagnostics(diag, "dummy", "dummy.h5");

    Diagnostics_loss d2();
    //b.set_diag_loss_aperture(d2);

    b.checkout_particles();
    auto parts = b.get_local_particles();
    parts(122, 3) = 122.3;
    parts(23, 5) = 23.5;
    b.checkin_particles();

    ar(sim);
}

void load()
{
    std::ifstream is(ar_name());
    cereal::JSONInputArchive ar(is);

    auto sim = Bunch_simulator::create_empty_bunch_simulator();

    ar(sim);

    Bunch& b = sim.get_bunch();
    b.checkout_particles();
    auto parts = b.get_local_particles();

    std::cout << "122.3 = " << parts(122, 3) << "\n";
    std::cout << "23.5 = " << parts(23, 5) << "\n";
    std::cout << "1000.6 = " << parts(1000, 6) << "\n";
    std::cout << "1024.2 = " << parts(1024, 2) << "\n";
}




int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    //run();

    save();
    load();


    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

