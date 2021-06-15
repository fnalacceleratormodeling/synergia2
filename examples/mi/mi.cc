
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/pcg_distribution.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/checkpoint.h"
#include "synergia/simulation/lattice_simulator.h"

#include "synergia/simulation/populate_stationary.h"
#include "synergia/simulation/split_operator_stepper.h"
//#include "synergia/simulation/independent_stepper_elements.h"

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
//#include "synergia/collective/space_charge_3d_open_hockney.h"
//#include "synergia/collective/space_charge_rectangular.h"
//#include "synergia/collective/dummy_collective_operator.h"

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

Lattice get_lattice()
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

void run()
{
    namespace LS = Lattice_simulator;

    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_STEP);

    auto lattice = get_lattice();

    LS::CourantSnyderLatticeFunctions(lattice);
    LS::calc_dispersions(lattice);

    for(auto const& elm : lattice.get_elements())
    {
        std::cout 
            << std::setprecision(15)
#if 0
            << elm.lf.beta.hor << ", " 
            << elm.lf.beta.ver << ", "
            << elm.lf.psi.hor << ", "
            << elm.lf.psi.ver << ", "
#endif
            << elm.lf.dispersion.hor << ", "
            << elm.lf.dispersion.ver << ", "
            << elm.lf.dPrime.hor << ", "
            << elm.lf.dPrime.ver << ", "

            << elm.lf.arcLength << ", "
            << "\n";
    }
}

void run_and_save(std::string & prop_str, std::string & sim_str)
{
    namespace LS = Lattice_simulator;

    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_STEP);


    // from opts
    double rf_voltage = 1.0/18; // RF cavity voltage in MV



    auto lsexpr = read_lsexpr_file("mi20_raw.lsx");
    Lattice lattice(lsexpr);

    lattice.set_all_string_attribute("extractor_type", "libff");

    screen << "lattice # elements: " << lattice.get_elements().size() << "\n";
    screen << "lattice length: " << lattice.get_length() << "\n";

    int harmno = 588;

    // closed orbit tolerance set to 1e-10
    LS::set_closed_orbit_tolerance(1e-6);

    // tune the lattice
    LS::tune_circular_lattice(lattice);

    // get the reference particle
    auto const & ref = lattice.get_reference_particle();
    double energy = ref.get_total_energy();
    double beta = ref.get_beta();
    double gamma = ref.get_gamma();

    screen(LV::INFO) 
        << "reference momentum = " << ref.get_momentum() << " GeV\n"
        << "energy: " << energy << "\n"
        << "beta: " << beta << "\n"
        << "gamma: " << gamma << "\n";

    // set rf cavity frequency
    // harmno * beta * c/ring_length
    double freq = harmno * beta * pconstants::c / lattice.get_length();

    screen 
        << "RF frequency: " << freq << "\n"
        << "Begin setting RF voltage... \n";

    // rf cavity voltage, is 1.0 MV total distributed over 18 cavities.  MAD8
    // expects cavities voltages in  units of MV.
    for (auto& elm : lattice.get_elements())
    {
        if (elm.get_type() == element_type::rfcavity)
        {
            elm.set_double_attribute("volt", rf_voltage);

            // set the harmonic number so the frequency is set 
            elm.set_double_attribute("harmon", harmno);

            // set the first pass frequency so I can get the bucket length
            elm.set_double_attribute("freq", freq*1.0e-6);
        }
    }

    screen << "Finish setting RF voltage...\n";

    auto tunes = LS::calculate_tune_and_cdt(lattice);
    auto chromes = LS::get_chromaticities(lattice);

    screen
        << "Unadjusted x tune: " << tunes[0] << "\n"
        << "Unadjusted y tune: " << tunes[1] << "\n"
        << "Unadjusted x chromaticity: " 
        << chromes.horizontal_chromaticity << "\n"
        << "Unadjusted y chromaticity: " 
        << chromes.vertical_chromaticity << "\n"
        ;


    // mark the focusing and defocusing elements for H/V tune correctors
    for (auto& elm : lattice.get_elements())
    {
        if (elm.get_type() == element_type::quadrupole)
        {
            auto ename = elm.get_name();

            if ( (ename.find("iqb") == 0) 
                    || (ename.find("iqc") == 0)
                    || (ename.find("iqd") == 0)
                    || (ename.find("iqe") == 0)
                    || (ename.find("iqd") == 0)
                    || (ename.find("iqg") == 0) )
            {
                double str = elm.get_double_attribute("k1", 0.0);

                if (str>0.0) elm.set_marker(marker_type::h_tunes_corrector);
                else if (str<0.0) elm.set_marker(marker_type::v_tunes_corrector);
            }
        }

        if (elm.get_type() == element_type::sextupole)
        {
            double str = elm.get_double_attribute("k2", 0.0);

            if (str>0.0) elm.set_marker(marker_type::h_chrom_corrector);
            else if (str<0.0) elm.set_marker(marker_type::v_chrom_corrector);
        }
    }

    // adjust tunes
    double xtune_adjust = 0.1;
    double ytune_adjust = 0.15;

    double xtune = tunes[0] + xtune_adjust;
    double ytune = tunes[1] + ytune_adjust;

    double tune_tolerance = 1.0e-8;

    screen
        << "Adjusting tunes to: \n"
        << "    xtune: " << xtune << "\n"
        << "    ytune: " << ytune << "\n"
        ;

    LS::adjust_tunes(lattice, xtune, ytune, tune_tolerance);

    // adjust chromaticities
    double xchrom_adjust = 0.1;
    double ychrom_adjust = 0.15;

    double xchrom = chromes.horizontal_chromaticity + xchrom_adjust;
    double ychrom = chromes.vertical_chromaticity + ychrom_adjust;

    screen
        << "Adjusting chromaticities to: \n"
        << "    xchrom: " << xchrom << "\n"
        << "    ychrom: " << ychrom << "\n"
        ;

    LS::adjust_chromaticities(lattice, xchrom, ychrom);

    // tunes and chromaticities after adjustments
    tunes = LS::calculate_tune_and_cdt(lattice);
    chromes = LS::get_chromaticities(lattice);

    screen
        << "Adjusted x tune: " << tunes[0] << "\n"
        << "Adjusted y tune: " << tunes[1] << "\n"
        << "Adjusted x chromaticity: " 
        << chromes.horizontal_chromaticity << "\n"
        << "Adjusted y chromaticity: " 
        << chromes.vertical_chromaticity << "\n"
        ;

     
    return;

    // space charge
    Space_charge_2d_open_hockney_options sc_ops(128, 128, 256);
    //Space_charge_rectangular_options sc_ops({64, 64, 128}, {0.1, 0.1, 1.0});
    sc_ops.comm_group_size = 1;

    // stepper
    Split_operator_stepper stepper(sc_ops, 71);

    // Propagator
    Propagator propagator(lattice, stepper);

    // bunch simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            //lattice.get_reference_particle(), 1024 * 1024 * 4, 2.94e10,
            //lattice.get_reference_particle(), 4194394, 2.94e10,
            lattice.get_reference_particle(), 1024*1024, 2.94e10,
            //lattice.get_reference_particle(), 16777216, 2.94e10,
            //lattice.get_reference_particle(), 0, 2.94e10,
            Commxx() );

    // get bunch
    auto & bunch = sim.get_bunch();

    // reserve particle slots
    //bunch.reserve(6000000);


    double emitx = 12.57e-6;
    double emity = 9.3e-6;
    double dpop = 2.4e-4/3.0;

    auto map = LS::get_linear_one_turn_map(lattice);

    auto map_x = Kokkos::subview(map, std::make_pair(0, 2), std::make_pair(0, 2));
    auto rx = LS::map_to_twiss(map_x);

    auto map_y = Kokkos::subview(map, std::make_pair(2, 4), std::make_pair(2, 4));
    auto ry = LS::map_to_twiss(map_y);

    auto map_z = Kokkos::subview(map, std::make_pair(4, 6), std::make_pair(4, 6));
    auto rz = LS::map_to_twiss(map_z);

    screen 
        << "Lattice parameters\n"
        << "alpha_x: " << rx[0] << ", alpha_y: " << ry[0] << "\n"
        << "beta_x:  " << rx[1] << ", beta_y:  " << ry[1] << "\n"
        << "q_x: " << rx[2] << ", q_y: " << ry[2] << ", q_s: " << rz[2] << "\n"
        << "beta_cdt: " << rz[1] << ", beta_longitudinal: " << rz[1]*beta << "\n";


    double bx = rx[1];
    double by = ry[1];
    double b_cdt = rz[1];

    double stdx = sqrt(bx*emitx/4.0);
    double stdy = sqrt(by*emity/4.0);
    double std_cdt = dpop*b_cdt;

    double half_bucket_length = LS::get_bucket_length(lattice) / 2.0;

    screen
        << "dpop = " << dpop << ", std_cdt = " << std_cdt << "\n\n";

    screen
        << "Beam parameters for bunch generation\n"
        << "invariant emittance_H: " << emitx << " [m-rad]\n"
        << "invariant emittance_V: " << emity << " [m-rad]\n"
        << "Delta-p/p RMS: " << dpop << "\n"
        << "RMS x: " << stdx << " [m]\n"
        << "RMS y: " << stdy << " [m]\n"
        << "RMS z: " << std_cdt*beta << " [m]"
        << ", RMS c*dt: " << std_cdt << " [m]"
        << ", RMS dt: " << std_cdt*1e9/pconstants::c << " [ns]\n"
        << "half_bucket_length: " << half_bucket_length << " [m]\n";

    // normal form for the lattice
    auto nf = LS::calculate_normal_form<2>(lattice);

    // actions
    auto actions = nf.stationaryActions(stdx, stdy, std_cdt);

    // min-max cdt
    double min_cdt = -half_bucket_length / beta;
    double max_cdt =  half_bucket_length / beta;

    // populate
    PCG_random_distribution dist(5);
    populate_6d_stationary_clipped_longitudinal_gaussian(
            dist, bunch, actions, min_cdt, max_cdt, nf);

#if 0
    // populate particle data
    karray1d means("means", 6);
    for (int i=0; i<6; ++i) means(i) = 0.0;

    karray2d_row covariances("covariances", 6, 6);
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


#if 0
    // or read from file
    //bunch.read_file_legacy("turn_particles_0000.h5");
    //bunch.write_file("bunch_particles_4M.h5");
    bunch.read_file("bunch_particles_4M.h5");

    //bunch.read_file_legacy("turn_particles16M_0000.h5");
    //bunch.write_file("bunch_particles_16M.h5");
    //bunch.read_file("bunch_particles_16M.h5");
#endif

    // statistics before propagate
    print_statistics(bunch, screen);

    Diagnostics_full2 diag_full2("diag_full2.h5");
    auto diag = sim.reg_diag_per_turn(diag_full2);

    Diagnostics_bulk_track diag_bt("diag_bulk_track.h5", 1000, 0);
    sim.reg_diag_per_turn(diag_bt);

    Diagnostics_particles diag_part("diag_particles.h5", 100);
    sim.reg_diag_per_turn(diag_part);

    sim.reg_diag_loss_aperture("loss.h5");

    // propagate
    propagator.propagate(sim, simlog, 2);

    // statistics after propagate
    print_statistics(bunch, screen);
    simple_timer_print(screen);

    // dump
    prop_str = propagator.dump();
    sim_str = sim.dump();

    //propagator.print_steps(screen);
    syn::checkpoint_save(propagator, sim);

    return;
}

void resume_and_save(std::string & prop_str, std::string & sim_str)
{
    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_STEP);

    auto propagator = Propagator::load_from_string(prop_str);
    auto sim = Bunch_simulator::load_from_string(sim_str);

    auto & bunch = sim.get_bunch();
    print_statistics(bunch, screen);

    //propagator.print_steps(screen);
    propagator.propagate(sim, simlog, 1);

    // timer and statistics
    simple_timer_print(screen);

    prop_str = propagator.dump();
    sim_str = sim.dump();
}

void checkpoint_resume()
{
    Logger screen(0, LV::DEBUG);
    Logger simlog(0, LV::INFO_STEP);

    screen << "resuming from checkpoint...\n";

    auto cp = syn::checkpoint_load();

    Propagator& propagator = cp.first;
    Bunch_simulator& sim = cp.second;

    auto & bunch = sim.get_bunch();
    print_statistics(bunch, screen);

    //propagator.print_steps(screen);
    propagator.propagate(sim, simlog, 2);

    // timer and statistics
    simple_timer_print(screen);
}


std::string ar_name()
{
    std::stringstream ss;
    ss << "cp-" << Commxx().rank() << ".json";
    return ss.str();
}

void bs_save()
{
    std::ofstream os(ar_name());
    cereal::JSONOutputArchive ar(os);

    auto sim = Bunch_simulator::create_single_bunch_simulator(
            Reference_particle(), 4194394, 2.94e10,
            Commxx() );

    Bunch& b = sim.get_bunch();

    Diagnostics_dummy diag;
    b.add_diagnostics(diag);

    Diagnostics_loss d2;
    //b.set_diag_loss_aperture(d2);

    b.checkout_particles();
    auto parts = b.get_local_particles();
    parts(122, 3) = 122.3;
    parts(23, 5) = 23.5;
    b.checkin_particles();

    ar(sim);
}

void bs_load()
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

//#include "synergia/utils/json.h"

int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    //run();
    //return 0;

    std::string prop_str, sim_str;
    run_and_save(prop_str, sim_str);

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

