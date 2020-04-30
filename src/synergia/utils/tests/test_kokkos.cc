#include "synergia/utils/catch.hpp"
#include "Kokkos_Core.hpp"

TEST_CASE("kokkos_padding", "[BunchParticles]")
{
    auto alloc = Kokkos::view_alloc("test", Kokkos::AllowPadding);
    auto darr = Kokkos::View<double**, Kokkos::LayoutLeft>(alloc, 11, 7);
    auto harr = Kokkos::create_mirror_view(darr);

    REQUIRE(darr.stride(0) == harr.stride(0));
    REQUIRE(darr.stride(1) == harr.stride(1));

    CHECK_NOTHROW(Kokkos::deep_copy(harr, darr));
}

