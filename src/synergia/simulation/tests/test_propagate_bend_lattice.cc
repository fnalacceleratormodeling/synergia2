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

cell: sequence, l=pi/2, refer=entry;
bpiover2, at=0.0;
endsequence;

circle: sequence, l=2*pi, refer=entry;
cell, at=0.0;
cell, at=pi/2;
cell, at=pi;
cell, at=3*pi/2;
endsequence;
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("circle");
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
    
	double cTref = 4 * L1/beta;
#if 0
	std::cout << "L1: " << std::setprecision(16) << L1 << std::endl;
	std::cout << "beta: " << std::setprecision(16) << beta <<std::endl;
	std::cout << "cTref: " << std::setprecision(16) << cTref << std::endl;
#endif

	// need beta for off-momentum particle
	double betagamma2 = beta * gamma * (1+dpp);
	double gamma2 = std::sqrt(betagamma2*betagamma2 + 1.0);
	double beta2 = betagamma2/gamma2;

	double cToff_momentum = 4 * L2/beta2;

#if 0
	std::cout << "L2: " << std::setprecision(16) << L2 << std::endl;
	std::cout << "beta2: " << std::setprecision(16) << beta2 <<std::endl;
	std::cout << "cToff_momentum: " << std::setprecision(16) << cToff_momentum << std::endl;
#endif

	double cTdiff = cToff_momentum - cTref;
	CHECK(bp(1, 0) == Approx(dpp));
	CHECK(cTdiff == Approx(bp(1, 4)));

	double cdt0 = Lattice_simulator::calculate_cdt(lattice, 0.0);
	CHECK( cdt0 == Approx(cTref) );

#if 0
	std::cout << std::setprecision(16) << "cdt0 from calculate_cdt" <<
		cdt0 << std::endl;
#endif

	double cdt1 = Lattice_simulator::calculate_cdt(lattice, dpp);
	CHECK( cdt1 == Approx(cToff_momentum) );

#if 0
	std::cout << std::setprecision(16) << "cdt1 from calculate_cdt" <<
		cdt1 << std::endl;
#endif

	auto chrom = Lattice_simulator::get_slip_factors(lattice);

	double slip_factor_geom = (cToff_momentum/cTref-1.0)/dpp;
	// chrom.slip_factor is calculated removing higher order corrections
	CHECK (chrom.slip_factor == Approx(slip_factor_geom).epsilon(4.0e-4) );

	std::cout << "slip factor from LS: " << std::setprecision(16) << chrom.slip_factor << std::endl;
	std::cout << "slip factor from geometry: " << std::setprecision(16) << slip_factor_geom  << std::endl;

}

