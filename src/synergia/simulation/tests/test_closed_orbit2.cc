#include "synergia/utils/catch.hpp"

#include "synergia/simulation/lattice_simulator.h"

#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/utils.h"
#include <iomanip>
#include <iostream>
#include <string>


Lattice
get_lattice()
{
    static std::string fodo_madx(R"foo(
ncells=4;
beam, particle=proton,energy=pmass+0.8;

bpiover2: sbend, l=pi/2, angle=pi/2;

cell: sequence, l=4+pi/2, refer=entry;
bpiover2, at=0.0;
endsequence;

square: sequence, l=16+2*pi, refer=entry;
cell, at=0.0;
cell, at=4.0+pi/2;
cell, at=8.0+pi;
cell, at=12.0+3*pi/2;
endsequence;
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("square");
}

static double cdt0 = 0.0;
static double cdt1 = 0.0;
static double beta0 = 0.0;
static double beta1 = 0.0;


TEST_CASE("closed_orbit_at_0dpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

#if 0
	std::cout << "lattice" << std::endl;
	std::cout << "<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	std::cout << lattice.as_string() << std::endl;
	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
#endif

    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice);

#if 0
    std::cout << "zero particle closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }
#endif

    for (int i=0; i<6; ++i) {
        CHECK (std::abs(closed_orbit_state[i]) < 1.0e-12);
    }

    // Check the c dT while we're at it.
    double cdt = Lattice_simulator::calculate_cdt(lattice);
#if 0
    std::cout << "on-momentum cdt: " << std::setprecision(16) << cdt << std::endl;
#endif

	cdt0 = cdt; // save on-momentum cdt for slip factor calculation

    double beta = lattice.get_reference_particle().get_beta();
	beta0 = beta;
    double length = lattice.get_length();
    CHECK( cdt == Approx(length/beta) );
    // calculate based on what I know the lattice is
    double L = 4*(4.0 + Kokkos::numbers::pi/2);
    CHECK( cdt == Approx(L/beta));

#if 0
    std::cout << "on-momentum geometric orbit length: " <<
             std::setprecision(16) << L << std::endl;
    std::cout << "on-momentum beta*calculate_cdt(): " << std::setprecision(16) <<  cdt*beta << std::endl;
#endif

}

TEST_CASE("closed_orbit_nonzerodpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    constexpr double dpp=1.0e-3;

    // determining the closed orbit state doesn't converge because this
    // lattice is unstable as defined, giving it a hint by setting
    // the reference particle state allows it to find one.

    std::array<double, 6> init_guess({dpp, 0.0, 0.0, 0.0, 0.0, dpp});

	// need to set initial closed orbit guess because this lattice is unstable
    lattice.get_reference_particle().set_state(init_guess);

	constexpr double co_tolerance = 1.0e-10;
	Lattice_simulator::set_closed_orbit_tolerance(co_tolerance);

    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice, dpp);

#if 0
    std::cout << "dpp=" << dpp << " off-mmentum closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }
#endif

    for (int i=0; i<4; ++i) {
        CHECK (std::abs(closed_orbit_state[i]-init_guess[i]) < co_tolerance);
    }

    // Check the c dT while we're at it.
    double cdt = Lattice_simulator::calculate_cdt(lattice, dpp);
	cdt1 = cdt; // save off-momentum cdt for slip factor

#if 0
    std::cout << "off-momentum  calculate_cdt: " << cdt << std::endl;
#endif

    // this is the off-momentum particle so beta is not the reference beta
    double p = lattice.get_reference_particle().get_momentum()*(1+dpp);
    double mass = lattice.get_reference_particle().get_mass();
    double betagamma = p/mass;
    double gamma = std::sqrt(betagamma*betagamma + 1);
    double beta = betagamma/gamma;
	beta1 = beta;

#if 0
    std::cout << "beta for off-momentum particle: " << beta << std::endl;
#endif
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

#if 0
	std::cout << "L1: " << L1 << std::endl;
	std::cout << "L2: " << L2 << std::endl;
	std::cout << "R1: " << R1 << std::endl;
	std::cout << "R2: " << R2 << std::endl;
#endif

    // 4 bends of length L2 plus 4 straights of lenght 4.0
	// This level of accuracy is not great
	CHECK( cdt == Approx(4*(4+L2)/beta).epsilon(4.0e-4) );

#if 0
    std::cout << std::setprecision(16) << "off-momentum geometric orbit length: " << 4*(4+L2) << std::endl;
    std::cout << std::setprecision(16) << "off-momentum beta* calculate_cdt(): " << cdt*beta << std::endl;
#endif

	//  get_slip_factors calls get_cdt() which calls calculate_closed_orbit
	// which fails off-momentum if the closed orbit state is not set to a good
	// guess. Does getting the 0dpp closed orbit work if it is set?

	lattice.get_reference_particle().set_state(
		{
			dpp, 0.0, 0.0, 0.0, 0.0, dpp
		}
		);
	auto zerodppcos = Lattice_simulator::calculate_closed_orbit(lattice, 0.0);

#if 0
	std::cout << "closed orbit state for 0dpp with initial off-momentum guess:" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << zerodppcos[i] << std::endl;    
    }
#endif

	double cdt_on_momentum_starting_with_guess = Lattice_simulator::calculate_cdt(lattice, 0.0);

#if 0
	std::cout << "cdt calculated starting from off-momentum closed orbit state: " << std::setprecision(16) << cdt_on_momentum_starting_with_guess << std::endl;
#endif

    // check cdt for the on-momentum closed orbit
	double cdt_on_momentum = 4 * (4 + L1)/beta0;
	CHECK( cdt_on_momentum_starting_with_guess == Approx(cdt_on_momentum).epsilon(4.0e-4) );

	double cdt_plus = Lattice_simulator::calculate_cdt(lattice, dpp);
	double cdt_minus = Lattice_simulator::calculate_cdt(lattice, -dpp);

#if 0
	std::cout << "cdt(+dpp) calculated starting from off-momentum closed orbit state" << std::setprecision(16) << cdt_plus << std::endl;
	std::cout << "cdt(-dpp) calculated starting from off-momentum closed orbit state" << std::setprecision(16) << cdt_minus << std::endl;
#endif

	// Check cdt for the off-momentum closed orbit
	double cdt_off_momentum = 4 * (4 + L2)/beta1;
	CHECK( cdt_plus == Approx(cdt_off_momentum).epsilon(4.0e-4) );

#if 0
	std::cout << "cdt calculated starting with 0 closed orbit state: " << std::setprecision(16) << cdt0 << std::endl;
#endif
	
	auto chrom = Lattice_simulator::get_slip_factors(lattice);

#if 0
	std::cout << "slip factor: " << chrom.slip_factor << std::endl;
#endif

	double slip_factor_should_be = ( (4 + L2)/(4 + L1) * (beta0/beta1) - 1.0)/dpp;
	
#if 0
	std::cout << "slip factor from LS: " << std::setprecision(16) << (1.0 - cdt1/cdt0)/dpp << std::endl;
	std::cout << "slip factor from geometry: " << std::setprecision(16) << slip_factor_should_be << std::endl;
	std::cout << "slip factor by hand from LS: " << std::setprecision(16) << (cdt_plus - cdt_minus)/(2*cdt_on_momentum_starting_with_guess*dpp) << std::endl;
#endif
	
	// chrom.slip_factor removes higher order contributions
	CHECK( slip_factor_should_be == Approx(chrom.slip_factor).epsilon(1.0e-2));
	

}
