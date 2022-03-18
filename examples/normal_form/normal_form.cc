
#include "synergia/lattice/madx_reader.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include <iostream>

const int order = 3; // higher orders use more memory and may not compile on many systems

void run()
{
    Logger screen(0, LoggerV::DEBUG);
    Logger log(0, LoggerV::INFO_TURN);

    // read lattice
    MadX_reader reader;
    auto lattice = reader.get_lattice("fodo", "channel.madx");

    // tune the lattice
    lattice.set_all_string_attribute("extractor_type", "libff");
    Lattice_simulator::tune_circular_lattice(lattice);

    lattice.print(screen);

    // calculate normal form
    auto nf = Lattice_simulator::calculate_normal_form<order>(lattice);

#if 0
    auto nfd = nf.cnvDataToNormalForm({0.1, 0.15, 0.2, 0.25, 0.05, 0.01});
    auto hfd = nf.cnvDataFromNormalForm(nfd);

    std::cout << "normal form = \n";
    for(int i=0; i<6; ++i) std::cout << nfd[i] << "\n";
    std::cout << "\n";

    std::cout << "human form = \n";
    for(int i=0; i<6; ++i) std::cout << hfd[i] << "\n";
    std::cout << "\n";

    return;
#endif

    // reference
    auto const & ref = lattice.get_reference_particle();

    auto energy = ref.get_total_energy();
    auto momentum = ref.get_momentum();
    auto gamma = ref.get_gamma();
    auto beta = ref.get_beta();

    int macro_particles = 80;
    double real_particles = 1.0e9;

    Independent_stepper_elements stepper(1);
    Propagator propagator(lattice, stepper);

    auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, macro_particles, real_particles, Commxx() );

    auto & bunch = sim.get_bunch();
    auto hparts = bunch.get_host_particles();

    // init particles
    for(int i=0; i<macro_particles; ++i)
        for(int j=0; j<6; ++j)
            hparts(i,j) = 0.0;

    for(int i=0; i<macro_particles/2; ++i)
        hparts(i, 0) = 0.001 * i;

    for(int i=macro_particles/2; i<macro_particles; ++i)
        hparts(i, 2) = 0.001 * (i - macro_particles/2);

    // check in
    bunch.checkin_particles();

    // nf data
    std::ofstream of("nf.dat");

    // propagate actions
    sim.reg_prop_action_turn_end([&screen, &nf, &of](
            auto& sim, auto& lattice, int turn) {
        //screen << "turn = " << turn << "\n";
        auto& bunch = sim.get_bunch();
        auto hparts = bunch.get_host_particles();
        auto nparts = bunch.size();

        bunch.checkout_particles();

        for(int i=0; i<nparts; ++i)
        {
            std::array<double, 6> coord;
            for(int j=0; j<6; ++j) coord[j] = hparts(i,j);

            auto nfd = nf.cnvDataToNormalForm(coord);

            of << nfd[0].real() << " "
               << nfd[0].imag() << " "
               << nfd[1].real() << " " 
               << nfd[1].imag() << " " 
               << nfd[2].real() << " " 
               << nfd[2].imag() << " ";
        }

        of << "\n";
        of.flush();
    });

    propagator.propagate(sim, log, 1000);

    bunch.checkout_particles();
    bunch.print_particle(0, screen);
    bunch.print_particle(1, screen);
    bunch.print_particle(2, screen);
    bunch.print_particle(3, screen);
    bunch.print_particle(40, screen);
    bunch.print_particle(41, screen);
    bunch.print_particle(42, screen);
    bunch.print_particle(43, screen);
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

