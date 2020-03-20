#include "synergia/utils/catch.hpp"

#include "synergia/bunch/bunch_particles.h"


TEST_CASE("BunchParticles", "[BunchParticles]")
{
    CHECK(true);

    BunchParticles bp("test", 11, -1, Commxx::World);
    std::cout << bp.size() << ", " << bp.capacity() << "\n";

    bp.reserve(17, Commxx::World);
    std::cout << bp.size() << ", " << bp.capacity() << "\n";
}
