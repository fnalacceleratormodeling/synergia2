
#include "synergia/lattice/madx_reader.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include <iostream>


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

    // one turn map
    auto map = Lattice_simulator::get_one_turn_map(lattice);
    std::cout << map.to_json() << "\n";

#if 0
    // reference
    auto const & ref = lattice.get_reference_particle();

    auto energy = ref.get_total_energy();
    auto momentum = ref.get_momentum();
    auto gamma = ref.get_gamma();
    auto beta = ref.get_beta();
#endif
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

