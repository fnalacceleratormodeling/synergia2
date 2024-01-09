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


// Put in focussing to stabilize lattice
Lattice
get_lattice()
{
    static std::string fodo_madx(R"foo(
ncells=4;
beam, particle=proton,energy=pmass+0.8;

bpiover2: sbend, l=pi/2, angle=pi/2;
fo: multipole, knl={0, 0.75};
defo: multipole, knl={0, -0.75};

cell: sequence, l=1+pi/2, refer=entry;
bpiover2, at=0.0;
fo, at=1/3;
defo, at=2/3;
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

TEST_CASE("closed_orbit_at_0dpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice);
    std::cout << "zero particle closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }

    for (int i=0; i<6; ++i) {
        CHECK (std::abs(closed_orbit_state[i]) < 1.0e-12);
    }
}

TEST_CASE("closed_orbit_nonzerodpp")
{
    Logger screen(0, LoggerV::INFO);

    Lattice lattice = get_lattice();

    constexpr double dpp=1.0e-3;

    auto closed_orbit_state = Lattice_simulator::calculate_closed_orbit(lattice, dpp);
    std::cout << "dpp=" << dpp << " off-mmentum closed orbit state" << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << std::setprecision(17) << i << ": " << closed_orbit_state[i] << std::endl;    
    }

}
