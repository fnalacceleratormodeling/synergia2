#define CATCH_CONFIG_RUNNER
#include "synergia/utils/catch.hpp"

#include "synergia/simulation/bunch_simulator.h"

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 100;
const double real_num = 2.0e12;
auto const id = Bunch::id;

TEST_CASE("Bunch_simulator", "[Bunch_simulator]")
{
    CHECK(true);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    int result = Catch::Session().run(argc, argv);

    Kokkos::finalize();
    MPI_Finalize();
    return result;
}
