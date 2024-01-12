#include "synergia/utils/catch.hpp"

#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/bunch/bunch.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/utils.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include <iomanip>
#include <iostream>
#include <string>


using LV = LoggerV;


Lattice
get_lattice()
{
    static std::string fodo_madx(R"foo(
ncells=4;
beam, particle=proton,energy=pmass+0.8;

bpiover2: sbend, l=pi/2, angle=pi/2;

cell: sequence, l=1+pi/2, refer=entry;
bpiover2, at=0.0;
endsequence;

square: sequence, l=4+2*pi, refer=entry;
cell, at=0.0;
cell, at=1.0+pi/2;
cell, at=2.0+pi;
cell, at=3.0+3*pi/2;
endsequence;
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("square");
}

Propagator get_propagator(Lattice const& lattice)
{
    Independent_stepper_elements stepper(1);
    Propagator prop(lattice, stepper);
    return prop;
}

TEST_CASE("closed_orbit_at_0dpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    auto refpart = lattice.get_reference_particle();

    Propagator prop(get_propagator(lattice));

    auto sim = Bunch_simulator::
        create_single_bunch_simulator(lattice.get_reference_particle(),
                                     16, 5.0e10, Commxx());
    
    // setup the bunch
    auto& bunch = sim.get_bunch();
    
    double beta = refpart.get_beta();
    double gamma = refpart.get_gamma();

    // calculate the closed orbit for on and off-momentum particles to
    // propagate

    /*
    # L = R theta => R = L/theta
    # R is radius of curvature of bend magnet ~ p/eB
    # R2/R1 = p2/p1
    R1 = bend.get_length()/bend.get_bend_angle()
    R2 = R1 * (1+dpp)

    L1 = R1 * np.pi/2
    L2 = R2 * np.arccos(1 - R1/R2)
    */

    constexpr double dpp = 1.0e-3; // offset momentum

    double R1 = 1.0;
    double R2 = R1 * (1+dpp);
    double L1 = R1 * Kokkos::numbers::pi/2; // pathlength in bend on-momentum
    double L2 = R2 * Kokkos::numbers::pi/2; // pathlength in bend off-momentum
    
    bunch.checkout_particles();

    auto bp = bunch.get_local_particles();

    for (int i=0; i<16; ++i) {
        for (int j=0; j<6; ++j) {
            bp(i, j) = 0.0;
        }
    }

    // particle 0 in on-momentum down the center
    // particle 1 is offset momentum with x transverse offset to match
    bp(1, 0) = dpp;
    bp(1, 5) = dpp;

    bunch.checkin_particles();

    Logger simlog(0, LV::INFO_STEP);

    prop.propagate(sim, simlog, 1); // single pass

    bunch.checkout_particles();

#if 0
    // print particles

    for (int i=0; i<2; ++i) {
        for (int j=0; j<6; ++j) {
            if (j != 0) std::cout << ", ";
            std::cout << std::setprecision(16) << bp(i, j);
        }
        std::cout << "\n" << std::endl;
    }                
#endif

    // check the particle coordinates
    // particle 0 should remain at 0
    for (int i=0; i<6; ++i) {
        CHECK(abs(bp(0, i)) < 1.0e-12);
    }

    // check the cdt of particle 1
    

}



#if 0
    std::cout << "zero particle closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }

    for (int i=0; i<6; ++i) {
        CHECK (std::abs(closed_orbit_state[i]) < 1.0e-12);
    }

    // Check the c dT while we're at it.
    double cdt = Lattice_simulator::calculate_cdt(lattice);
    std::cout << "on-momentum cdt: " << std::setprecision(16) << cdt << std::endl;
    double beta = lattice.get_reference_particle().get_beta();
    double length = lattice.get_length();
    CHECK( cdt == Approx(length/beta) );
    // calculate based on what I think the lattice is
    double L = 4*(1 + Kokkos::numbers::pi/2);
    CHECK( cdt == Approx(L/beta));

    std::cout << "on-momentum geometric orbit length: " <<
             std::setprecision(16) << L << std::endl;
    std::cout << "on-momentum beta*calculate_cdt(): " << std::setprecision(16) <<  cdt*beta << std::endl;

}
#endif

#if 0
TEST_CASE("closed_orbit_nonzerodpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    constexpr double dpp=1.0e-3;

    // determining the closed orbit state doesn't converge because this
    // lattice is unstable as defined, but we can give it a hint by setting
    // the reference particle state.

    std::array<double, 6> init_guess({dpp, 0.0, 0.0, 0.0, 0.0, dpp});

    lattice.get_reference_particle().set_state(init_guess);

    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice, dpp);
    std::cout << "dpp=" << dpp << " off-mmentum closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }

    for (int i=0; i<4; ++i) {
        CHECK (std::abs(closed_orbit_state[i]-init_guess[i]) < 1.0e-12);
    }

    // Check the c dT while we're at it.
    double cdt = Lattice_simulator::calculate_cdt(lattice, dpp);
    std::cout << "off-momentum  calculate_cdt: " << cdt << std::endl;

    // this is the off-momentum particle so beta is not the reference beta
    double p = lattice.get_reference_particle().get_momentum()*(1+dpp);
    double mass = lattice.get_reference_particle().get_mass();
    double betagamma = p/mass;
    double gamma = std::sqrt(betagamma*betagamma + 1);
    double beta = betagamma/gamma;
    std::cout << "beta for off-momentum particle: " << beta << std::endl;
    double length = lattice.get_length();
    // calculate expanded pathlength with dpp

    /*
    # L = R theta => R = L/theta
    # R is radius of curvature of bend magnet ~ p/eB
    # R2/R1 = p2/p1
    R1 = bend.get_length()/bend.get_bend_angle()
    R2 = R1 * (1+dpp)

    L1 = R1 * np.pi/2
    L2 = R2 * np.arccos(1 - R1/R2)
    */

    double R1 = 1.0;
    double R2 = R1 * (1+dpp);
    double L1 = R1 * Kokkos::numbers::pi/2; // pathlength in bend on-momentum
    double L2 = R2 * Kokkos::numbers::pi/2; // pathlength in bend off-momentum

std::cout << "L1: " << L1 << std::endl;
std::cout << "L2: " << L2 << std::endl;
std::cout << "R1: " << R1 << std::endl;
std::cout << "R2: " << R2 << std::endl;

    // 4 bends of length L2 plus 4 straights of lenght 1.0
    CHECK( cdt == Approx(4*(1+L2)/beta) );
    std::cout << std::setprecision(16) << "off-momentum geometric orbit length: " << 4*(1+L2) << std::endl;
    std::cout << std::setprecision(16) << "off-momentum beta* calculate_cdt(): " << cdt*beta << std::endl;
}

#endif
