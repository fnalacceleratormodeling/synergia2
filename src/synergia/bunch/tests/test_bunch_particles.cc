#include "synergia/utils/catch.hpp"

#include "synergia/bunch/bunch_particles.h"

// number of particles
const int np = 11;

// idx of lost particles
const auto losts = {3, 5};

void init_particle_values(BunchParticles& bp)
{
    for(int i=0; i<np; ++i)
        for(int j=0; j<6; ++j)
            bp.hparts(i, j) = i+j*0.1;

    for(auto l : losts) bp.hmasks(l) = 0;

    bp.checkin_particles();
    bp.update_valid_num();
}

void check_particle_values(BunchParticles const& bp)
{
    bp.checkout_particles();

    for(int i=0; i<np; ++i)
    {
        for(int j=0; j<6; ++j)
            CHECK( bp.hparts(i, j) == i+j*0.1 );

        auto it = std::find(losts.begin(), losts.end(), i);
        if (it == losts.end()) CHECK(bp.hmasks(i) == 1);
        else                   CHECK(bp.hmasks(i) == 0);
    }
}

TEST_CASE("BunchParticles", "[BunchParticles]")
{
    // create object
    BunchParticles bp("test", np, -1, Commxx::World);

    // init particle data
    init_particle_values(bp);

    // check
    REQUIRE(bp.size() == np);
    REQUIRE(bp.capacity() >= np);
    REQUIRE(bp.num_valid() == np - losts.size());

    SECTION("verify values")
    {
        check_particle_values(bp);
    }

    SECTION("resizing bigger")
    {
        bp.reserve(np+6, Commxx::World);

        REQUIRE(bp.size() == np);
        REQUIRE(bp.capacity() >= np+6);
        REQUIRE(bp.num_valid() == np - losts.size());

        check_particle_values(bp);
    }

    SECTION("resizing smaller")
    {
        bp.reserve(np-6, Commxx::World);

        REQUIRE(bp.size() == np);
        REQUIRE(bp.capacity() >= np);
        REQUIRE(bp.num_valid() == np - losts.size());

        check_particle_values(bp);
    }
}
