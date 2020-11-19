

#include "synergia/lattice/lattice.h"

namespace Lattice_simulator
{
    constexpr const double default_closed_orbit_tolerance = 1.0e-13;

    // Both tune_linear_lattice() and tune_circular_lattice() set the frequency of the 
    // rfcavities based on the momentum of the lattice reference particle.  
    //
    // tune_linear_lattice uses the state of the reference particle as the starting 
    // point for propagation.  
    //
    // tune_circular_lattice() calculates a closed orbit and uses that as the starting 
    // point.  
    //
    // Both returns return an array of state, the transverse coordinates are the final 
    // coordinates. cdt is the c * the total propagation time of the particle which 
    // gives the total path length.

    // find the closed orbit for the circular lattice, propagate the lattice reference 
    // particle through the lattice slices using the closed orbit, and set the reference 
    // c*t for each lattice slice after tuning.
    // return values is the state for calcualted closed orbit
    // note that all the rf cavities will be set to 0 strength during the tuning process
    std::array<double, 6> 
    tune_linear_lattice(Lattice & lattice);

    std::array<double, 6> 
    tune_circular_lattice(Lattice & lattice, double tolerance=1.0e-13);

    std::array<double, 6> 
    tune_rfcavities(Lattice & lattice);

    // closed orbit
    std::array<double, 6> 
    calculate_closed_orbit(Lattice const& lattice, double dpp = 0.0, 
            double tolerance = default_closed_orbit_tolerance);

}
