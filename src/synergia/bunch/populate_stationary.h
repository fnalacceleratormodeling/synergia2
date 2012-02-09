#ifndef POPULATE_STATIONARY_H_
#define POPULATE_STATIONARY_H_

#include "synergia/foundation/distribution.h"
#include "synergia/bunch/bunch.h"
#include "synergia/simulation/lattice_simulator.h"


/// Populate a bunch with a shell of particles having fixed constant actions
/// but uniform in phase angle in all three planes in normal form space.
/// Any distribution or statistic constructed from phase space variables
/// will be stationary under propagation through the lattice (up to
/// limited statistics) if there is no other physics.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param I0 Action variable for 0 coordinate
/// @param I1 Action variable for 1 coordinate
/// @param I2 Action variable for 2 coordinate
void
populate_6d_stationary_torus(Distribution &dist, Bunch &bunch, double I0,
			     double I1, double I2, Lattice_simulator& lattice_simulator);


#endif /* POPULATE_STATIONARY_H_ */
