#include "Kokkos_Core_fwd.hpp"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/bunch_particles.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/catch.hpp"

#include <Kokkos_NumericTraits.hpp>
#include <Kokkos_Random.hpp>

constexpr double mass = 100.0;
constexpr double total_energy = 125.0;
constexpr int total_num = 100;
constexpr double real_num = 2.0e12;
constexpr auto id = Bunch::id;
constexpr int num_parts = 1024;

TEST_CASE("Bunch", "[Bunch]")
{
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);
    Bunch bunch(ref, num_parts, 1e13, Commxx());

    CHECK(bunch.get_local_num() == num_parts);

    bunch.checkout_particles();
    auto parts = bunch.get_host_particles();

    CHECK(parts(3, id) == 3);
    CHECK(parts(307, id) == 307);

    auto idx1 = bunch.search_particle(12);
    CHECK(idx1 == 12);

    auto idx2 = bunch.search_particle(13, 12);
    CHECK(idx2 == 13);

    auto idx3 = bunch.search_particle(15, 15);
    CHECK(idx3 == 15);

    auto p2 = bunch.get_particles_in_range_row(123, 6);
    CHECK(p2(0, 6) == 123);
    CHECK(p2(1, 6) == 124);
    CHECK(p2(4, 6) == 127);
}

#if defined SYNERGIA_HAVE_OPENPMD

void
check_particle_values(BunchParticles const& bp1, BunchParticles const& bp2)
{
    bp1.checkout_particles();

    for (int i = 0; i < num_parts; ++i) {
        for (int j = 0; j < 6; ++j) {
            CHECK(bp1.hparts(i, j) == bp2.hparts(i, j));
            CHECK(bp1.hmasks(i) == bp2.hmasks(i));
        }
    }
}

TEST_CASE("BunchI/O", "[Bunch]")
{
    Four_momentum fm(mass, total_energy);
    Reference_particle ref(pconstants::proton_charge, fm);
    Bunch bunch(ref, num_parts, 1e13, Commxx());

    CHECK(bunch.get_local_num() == num_parts);

    auto bp1 = bunch.get_bunch_particles();
    std::cout << "bp1 numactive is : " << bp1.num_active() << "\n";

    bp1.checkout_particles();
    bp1.reserve(num_parts, bunch.get_comm());

    Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>
        random_pool(/*seed=*/12345);
    Kokkos::fill_random(bp1.hparts, random_pool, 100.0);

    Kokkos::parallel_for(
        "test_bunch_fill_masks",
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, num_parts),
        KOKKOS_LAMBDA(const int& i) {
            // acquire the state of the random number generator engine
            auto generator = random_pool.get_state();
            double x = generator.drand(0., 100.);
            if (x > 50) {
                bp1.hmasks(i) = 1;
            } else {
                bp1.hmasks(i) = 0;
            }
            // do not forget to release the state of the engine
            random_pool.free_state(generator);
        });

    bp1.update_valid_num();
    bp1.update_total_num(bunch.get_comm());

    SECTION("write/read file")
    {
        {
            bunch.write_openpmd_file("test_bunch_dump.h5");
        }

        SECTION("read into 0 sized bunch particle")
        {
            Four_momentum fm(mass, total_energy);
            Reference_particle ref(pconstants::proton_charge, fm);
            Bunch bunch2(ref, 0, 0, Commxx());

            bunch2.read_openpmd_file("test_bunch_dump.h5");

            auto bp2 = bunch2.get_bunch_particles();
            bp2.checkout_particles();

            check_particle_values(bp1, bp2);
        }
    }
}
#endif
