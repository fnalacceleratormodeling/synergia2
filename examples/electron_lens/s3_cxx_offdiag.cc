
#include "synergia/foundation/physical_constants.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/checkpoint.h"
#include "synergia/simulation/lattice_simulator.h"

#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/bunch/populate.h"
#include "synergia/bunch/core_diagnostics.h"
//#include "synergia/bunch/diagnostics_loss.h"
//#include "synergia/bunch/diagnostics_track.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"

#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/lsexpr.h"

#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
//#include "synergia/simulation/dummy_collective_operator.h"

#include "synergia/utils/simple_timer.h"

using LV = LoggerV;

void print_statistics(Bunch & bunch, Logger & logger)
{
    logger(LV::DEBUG)
        << "Bunch statistics: "
        << "num_valid = " << bunch.get_local_num()
        << ", size = " << bunch.size()
        << ", capacity = " << bunch.capacity()
        << ", total_num = " << bunch.get_total_num()
        <<"\nMean and std: ";


    // print particles after propagate
    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std  = Core_diagnostics::calculate_std(bunch, mean);

    logger(LV::DEBUG) 
        << std::resetiosflags(std::ios::fixed)
        << std::setprecision(16)
        << std::setiosflags(std::ios::showpos | std::ios::scientific)
        << "\n"
        //<< "\nmean\tstd\n"
        ;

    for (int i=0; i<6; ++i) 
        logger(LV::DEBUG) << mean[i] << ", " << std[i] << "\n";

    logger(LV::DEBUG)
        << std::resetiosflags(std::ios::showpos | std::ios::scientific)
        << "\n";

    for (int p=0; p<4; ++p) bunch.print_particle(p, logger);

    logger << "\n";
}

void activate_elens(Lattice& lattice)
{
    const double elens_energy = 0.01; //opts.elensenergy;
    const double elens_radius = 0.0041525; //opts.elensradius;
    const double elens_length = 1.0; //opts.elenslength;
    const double elens_current = 0.0; //opts.elenscurrent;
    const double elens_longrms = 0.5; //opts.elenslongrms;

    const int elens_divide = 2; //opts.elensdivide;

    for(auto& ele : lattice.get_elements())
    {
        if (ele.get_type() != element_type::elens) continue;

        ele.set_string_attribute("extractor_type", "libff");
        ele.set_double_attribute("radius", elens_radius);
        ele.set_double_attribute("gaussian", 1.0);
        ele.set_double_attribute("eenergy", elens_energy);
        ele.set_double_attribute("l", 0.0);
        ele.set_double_attribute("current", elens_current*elens_length);
        ele.set_double_attribute("longrms", elens_longrms);
    }
}

void run_and_save(std::string & prop_str, std::string & sim_str)
{
    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_STEP);
    //Logger simlog(0, LV::DEBUG);

    auto lsexpr = read_lsexpr_file("adjusted_lattice.lsx");
    Lattice lattice(lsexpr);

    activate_elens(lattice);

    lattice.set_all_string_attribute("extractor_type", "libff");

    // tune the lattice
    //Lattice_simulator::tune_circular_lattice(lattice);

    // get the reference particle
    auto const & ref = lattice.get_reference_particle();

    screen(LV::INFO) 
        << "reference momentum = " << ref.get_momentum() << " GeV\n";

    // space charge
    Space_charge_3d_open_hockney_options sc_ops(64, 64, 128);
    sc_ops.green_fn = green_fn_t::linear;
    sc_ops.comm_group_size = 1;
    //sc_ops.periodic_z = false;
    //sc_ops.longitudinal_kicks = false;

    //Dummy_CO_options sc_ops;

    // stepper
    Split_operator_stepper stepper(72, sc_ops);
    //Split_operator_stepper stepper(sc_ops, 72);
    //stepper.append_collective_op(Dummy_CO_options());
    //Independent_stepper_elements stepper(1);

    // Propagator
    Propagator propagator(lattice, stepper);

    // bunch simulator
#if 0
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            //lattice.get_reference_particle(), 1024 * 1024 * 4, 2.94e10,
            //lattice.get_reference_particle(), 4194394, 2.94e10,
            //lattice.get_reference_particle(), 16777216, 2.94e10,
            lattice.get_reference_particle(), 0, 2.94e10,
            Commxx() );
#endif

    auto sim = Bunch_simulator::create_bunch_train_simulator(
            lattice.get_reference_particle(), 
            1000000, // opts.macroparticles, 
            2e11, // opts.real_particles,
            1, // opts.num_bunches, 
            5.759999999999998 // spacing, lattice_simulator->get_bucket_length()
            );

    // get bunch
    auto & bunch = sim.get_bunch();

    // TODO
    // for(bunch in bunches) {
    //   bunch.set_bucket_index(i)
    //   bunch.set_z_period_length(lattice.length() / opts.harmon);
    // }

    // reserve particle slots
    //bunch.reserve(6000000);

#if 1
    // populate particle data
    karray1d means("means", 6);
    for (int i=0; i<6; ++i) means(i) = 0.0;

    karray2d_row covariances("covariances", 6, 6);
    for (int i=0; i<6; ++i)
        for (int j=0; j<6; ++j)
            covariances(i, j) = 0.0;


    covariances(0,0) = 1.733566436318178e-05;
    covariances(0,1) = -1.851669292545274e-06;

    covariances(1,0) = -1.851669292545274e-06;
    covariances(1,1) = 2.555307588570028e-07;

    covariances(2,2) = 1.72342639987807e-05;
    covariances(2,3) = 1.903323834339802e-06;

    covariances(3,2) = 1.903323834339802e-06;
    covariances(3,3) = 2.682886740716003e-07;

    covariances(4,4) = 0.3527858128388799;
    covariances(4,5) = 3.351492383857787e-17;

    covariances(5,4) = 3.351492383857787e-17;
    covariances(5,5) = 8.29770226173816e-06;

    Random_distribution dist(13 /*opts.seed*/, bunch.get_comm());
    populate_6d(dist, bunch, means, covariances);
#endif


#if 0
    // or read from file
    //bunch.read_file_legacy("turn_particles_0000_4M.h5");
    //bunch.write_file("bunch_particles_4M.h5");
    bunch.read_file("bunch_particles_4M.h5");

    //bunch.read_file_legacy("turn_particles16M_0000.h5");
    //bunch.write_file("bunch_particles_16M.h5");
    //bunch.read_file("bunch_particles_16M.h5");
#endif

    // statistics before propagate
    print_statistics(bunch, screen);

#if 0
    Diagnostics_full2 diag_full2;
    sim.reg_diag_per_turn(diag_full2, "full2", "diag_full2.h5");

    Diagnostics_bulk_track diag_bt(1000, 0);
    sim.reg_diag_per_turn(diag_bt, "bulk_track", "diag_bulk_track.h5");

    Diagnostics_particles diag_part(100);
    sim.reg_diag_per_turn(diag_part, "particles", "diag_particles.h5");

    sim.reg_diag_loss_aperture("loss.h5");
#endif

    // propagate
    propagator.propagate(sim, simlog, 1);

    // statistics after propagate
    print_statistics(bunch, screen);
    simple_timer_print(screen);

    // dump
    prop_str = propagator.dump();
    sim_str = sim.dump();

    //propagator.print_steps(screen);
    //syn::checkpoint_save(propagator, sim);

    return;
}

int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    std::string prop_str, sim_str;
    run_and_save(prop_str, sim_str);

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

