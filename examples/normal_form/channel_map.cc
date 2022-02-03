
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
    auto mapping = Lattice_simulator::get_one_turn_map<3>(lattice);
    std::cout << mapping.to_json() << "\n";

#if 0
    // reference
    auto const & ref = lattice.get_reference_particle();

    auto energy = ref.get_total_energy();
    auto momentum = ref.get_momentum();
    auto gamma = ref.get_gamma();
    auto beta = ref.get_beta();
#endif

    for(int comp=0; comp<6; ++comp)
    {
        auto trigon = mapping[comp];
        std::cout << "\n\nComponent " << comp << "\n";

        trigon.each_term([&trigon](int i, auto const& inds, auto term) {
            std::cout << "power = " << inds.size();

            std::array<int, trigon.dim> exp{0};
            for(auto const& x : inds) ++exp[x];

            std::cout << ", exp = [";
            for(int x=0; x<exp.size()-1; ++x) 
                std::cout << exp[x] << ", ";
            std::cout << exp.back() << "]";

            std::cout << ", term = " << term << "\n";
        });
    }
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

