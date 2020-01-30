#include "synergia/utils/catch.hpp"

#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/bunch_simulator_impl.h"

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 100;
const double real_num = 2.0e12;
auto const id = Bunch::id;

const Reference_particle ref(1, 1.0, 1.0);

void check_simulator(Bunch_simulator const& sim, 
        int nb_pt, int nb_st, int size, int rank)
{
    CHECK(sim[0].get_num_bunches() == nb_pt);
    CHECK(sim[1].get_num_bunches() == nb_st);
}

TEST_CASE("create single bunch", "[Bunch_simulator]")
{
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            ref, 1024, 1e10, Commxx() );

    int mpi_size = Commxx().size();
    int mpi_rank = Commxx().rank();

    check_simulator(sim, 1, 0, mpi_size, mpi_rank);

    auto const & b = sim.get_bunch(0, 0);
}
