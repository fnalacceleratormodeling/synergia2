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
cell_len = 1;
beam, particle=proton,energy=pmass+0.8;

m1: marker;
m2: marker;
m3: marker;
m4: marker;

cell1: sequence, l=cell_len, refer=entry;
m1, at=cell_len;
endsequence;

cell2: sequence, l=cell_len, refer=entry;
m2, at=cell_len;
endsequence;

cell3: sequence, l=cell_len, refer=entry;
m3, at=cell_len;
endsequence;

cell4: sequence, l=cell_len, refer=entry;
m4, at=cell_len;
endsequence;

square: sequence, l=4*cell_len, refer=entry;
cell1, at=0.0;
cell2, at=cell_len;
cell3, at=2.0*cell_len;
cell4, at=3.0*cell_len;;
endsequence;
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("square");
}

TEST_CASE("closed_orbit_at_0dpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

#if 0
    std::cout << lattice.as_string() << std::endl;
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
    double beta = lattice.get_reference_particle().get_beta();
    double length = lattice.get_length();
    CHECK( cdt == Approx(length/beta) );
    // calculate based on what I think the lattice is
    double L = 4.0;
    CHECK( cdt == Approx(L/beta));

#if 0
    std::cout << "on-momentum orbit length: " <<
             std::setprecision(16) << L << std::endl;
#endif

}

TEST_CASE("closed_orbit_nonzerodpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();
    double L = lattice.get_length();

    constexpr double dpp=1.0e-3;

    // determining the closed orbit state doesn't converge because this
    // lattice is unstable as defined, but we can give it a hint by setting
    // the reference particle state.

    //std::array<double, 6> init_guess({dpp, 0.0, 0.0, 0.0, 0.0, dpp});

    //lattice.get_reference_particle().set_state(init_guess);

    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice, dpp);

#if 0
    std::cout << "dpp=" << dpp << " off-mmentum closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }
#endif

    for (int i=0; i<4; ++i) {
        CHECK (std::abs(closed_orbit_state[i]) < 1.0e-12);
    }

    // Check the c dT while we're at it.
    double cdt = Lattice_simulator::calculate_cdt(lattice, dpp);
    double p = lattice.get_reference_particle().get_momentum()*(1+dpp);
    double mass = lattice.get_reference_particle().get_mass();
    double betagamma = p/mass;
    double gamma = std::sqrt(betagamma*betagamma + 1);
    double beta = betagamma/gamma;

    CHECK( cdt == Approx(L/beta) );

#if 0
    std::cout << std::setprecision(16) << "off-momentum orbit length: " << L << std::endl;
    std::cout << std::setprecision(16) << "off-momentum length from cdt: " << cdt*beta << std::endl;
#endif

}
