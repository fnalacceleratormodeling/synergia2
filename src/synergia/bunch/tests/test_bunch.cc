#define CATCH_CONFIG_RUNNER
#include "synergia/utils/catch.hpp"

#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/bunch.h"

const double mass = 100.0;
const double total_energy = 125.0;
const int total_num = 100;
const double real_num = 2.0e12;


TEST_CASE("Bunch", "[Bunch]")
{
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);
    Bunch bunch(ref, 1024, 1e13, Commxx());

    CHECK(bunch.get_local_num() == 1024);

    bunch.checkout_particles();
    auto parts = bunch.get_host_particles();

    CHECK(parts(3, Bunch::id) == 3);
    CHECK(parts(307, Bunch::id) == 307);

    auto idx1 = bunch.search_particle(12);
    CHECK(idx1 == 12);

    auto idx2 = bunch.search_particle(13, 12);
    CHECK(idx2 == 13);

    auto idx3 = bunch.search_particle(15, 15);
    CHECK(idx3 == 15);
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
