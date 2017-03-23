#ifndef POPULATE_STATIONARY_H_
#define POPULATE_STATIONARY_H_

#include "synergia/foundation/distribution.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/fast_normal_form.h"


/// Populate a bunch with a shell of particles having fixed constant actions
/// but uniform in phase angle in all three planes in normal form space.
/// Any distribution or statistic constructed from phase space variables
/// will be stationary under propagation through the lattice (up to
/// limited statistics) if there is no other physics.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param actions std::vector<double> (3) the three mean actions
void
populate_6d_stationary_torus(Distribution &dist, Bunch &bunch, std::vector<double> actions, Lattice_simulator& lattice_simulator);

void
populate_6d_stationary_gaussian(Distribution &dist, Bunch &bunch, std::vector<double> actions, Lattice_simulator& lattice_simulator);

void
populate_6d_stationary_gaussian_adjust(Distribution &dist, Bunch &bunch, std::vector<double> actions, Lattice_simulator& lattice_simulator,
                                      Const_MArray1d_ref means, Const_MArray2d_ref covariances);
                                      

/// use limit_x =xmax/sigma_x  where xmax=aperture maximum  radius
void
populate_6d_stationary_gaussian_truncated(Distribution &dist, Bunch &bunch,
                const std::vector<double> actions,
                Lattice_simulator& lattice_simulator,   Const_MArray1d_ref limits );  

void
populate_6d_stationary_gaussian_truncated (Distribution &dist, Bunch &bunch,
                const std::vector<double> actions,
                Fast_normal_form & fnf ,   Const_MArray1d_ref limits); 

void
populate_6d_stationary_truncated_longitudinal_gaussian(Distribution &dist, Bunch &bunch, std::vector<double> actions, double n_sigma, Lattice_simulator& lattice_simulator);

void
populate_6d_stationary_clipped_longitudinal_gaussian(Distribution &dist,
                                                     Bunch &bunch,
                                                     const std::vector<double> actions,
                                                     double cdt_min, double cdt_max,
                                                     Lattice_simulator& lattice_simulator);


// helper functions for populate stationary*
void
fill_6d_normal_form_coords(Distribution &dist, MArray2d_ref nf_particles, const std::vector<double> actions);

void
fill_6d_normal_form_coords(Distribution &dist, MArray2d_view nf_particles, const std::vector<double> actions);

#endif /* POPULATE_STATIONARY_H_ */
