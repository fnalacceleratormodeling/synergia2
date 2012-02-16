#include "populate_stationary.h"

#include <iostream>

#include "synergia/foundation/math_constants.h"

void
populate_6d_stationary_torus(Distribution &dist, Bunch &bunch,
			     const std::vector<double> actions,
			     Lattice_simulator& lattice_simulator)
{
    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];

    // the normalized radii in normal form space
    double sr0 = sqrt(actions[0]);
    double sr1 = sqrt(actions[1]);
    double sr2 = sqrt(actions[2]);

    // fill even columns with random phases 0 - 2*pi
    for (int c=0; c<6; c+=2) {
      dist.fill_uniform(particles[boost::indices[range()][c]], 0.0, 2.0*mconstants::pi);
    }

    // go through and set the components in complex cartesian normal form space
    // eq. 13 of Theory and Praxis of Map Analysis
    // a_k = i \sqrt{I_k} e^{-i \phi_k }
    for (int part = 0; part < num_part; ++part) {
      double phase;
      phase = particles[part][0];
      particles[part][0] = sr0*sin(phase);
      particles[part][1] = -sr0*cos(phase);

      phase = particles[part][2];
      particles[part][2] = sr1*sin(phase);
      particles[part][3] = -sr1*cos(phase);

      phase = particles[part][4];
      particles[part][4] = sr2*sin(phase);
      particles[part][5] = -sr2*cos(phase);

    }

    lattice_simulator.convert_normal_to_human(particles);
}

void
populate_6d_stationary_gaussian(Distribution &dist, Bunch &bunch,
				const std::vector<double> actions,
				Lattice_simulator& lattice_simulator)
{
    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];

    // fill even columns with uniform [0,1) which will be converted
    // to exponential distribution later
    for (int c=0; c<6; c+=2) {
      dist.fill_uniform(particles[boost::indices[range()][c]], 0.0, 1.0);
    }

    // fill odd columns with random phases 0 - 2*pi
    for (int c=1; c<6; c+=2) {
      dist.fill_uniform(particles[boost::indices[range()][c]], 0.0, 2.0*mconstants::pi);
    }

    // go through and set the components in complex cartesian normal form space
    // eq. 13 of Theory and Praxis of Map Analysis
    // a_k = i \sqrt{I_k} e^{-i \phi_k }
    for (int c=0; c<3; ++c) {
      for (int part = 0; part < num_part; ++part) {
	// the action is drawn from an exponential distribution with mean
	// actions[c] according to -mean*log(1-r) where r is random number
	// [0,1).  We need the square root of the action for the normal form
	// coordinate according to Eq. 13.
	double square_root_action = sqrt(-actions[c]*log(1.0-particles[part][2*c]));
	double phase = particles[part][2*c+1];
	particles[part][2*c] = square_root_action * sin(phase);
	particles[part][2*c+1] = -square_root_action * cos(phase);
      }
    }

    lattice_simulator.convert_normal_to_human(particles);
}
