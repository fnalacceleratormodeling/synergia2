#include "populate_stationary.h"
#include "diagnostics.h"

#include <iostream>

#include "synergia/foundation/math_constants.h"

void
populate_6d_stationary_torus(Distribution &dist, Bunch &bunch,
			     const double I0,const double I1,
			     const double I2,
			     Lattice_simulator& lattice_simulator)
{
    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];

    std::cout << "in populate_6d_staionary_torus" << std::endl;
    std::cout << "bunch size: " << num_part << " x " <<
      particles.shape()[1] << std::endl;

    // the normalized radii in normal form space
    double sr0 = sqrt(I0);
    double sr1 = sqrt(I1);
    double sr2 = sqrt(I2);

    std::cout << "coordinate maximum extents: " <<
      "x: " << sr0 << "   y: " << sr1 << "  z: " << sr2 << std::endl;

    // fill even columns with random phases 0 - 2*pi
    for (int c=0; c<6; c+=2) {
      std::cout << "filling particles[:][" << c << "]" << std::endl;
#if 0
      dist.fill_uniform(particles[boost::indices[range()][c]], 0.0, 2.0*mconstants::pi);
#else
      double pstp = 2.0*mconstants::pi/num_part;
      double ph = 0.0;
      for (int i=0; i<num_part; ++i) {
	particles[i][c] = ph;
	ph += pstp;
      }
#endif
    }

    for (int c=0; c<6; c+=2) {
      double pmin=1.0e10;
      double pmax=-1.0e10;
      double psum=0.0;
      double psum2=0.0;
      for (int i=0; i<num_part; ++i) {
	if (particles[i][c] < pmin) {
	  pmin = particles[i][c];
	}
	if (particles[i][c] > pmax) {
	  pmax = particles[i][c];
	}
	psum += particles[i][c];
	psum2 += pow(particles[i][c],2);
      }
      std::cout << "statistics for particles[:][" << c << "]: " << std::endl;
      std::cout << "    min: " << pmin << std::endl;
      std::cout << "    max: " << pmax << std::endl;
      std::cout << "    mean: " << psum/particles.shape()[0] << std::endl;
      std::cout << "    std:  " << sqrt(psum2/particles.shape()[0] -
					pow( (psum/particles.shape()[0]), 2))
		<< std::endl;
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

    for (int part = 0; part < num_part; ++part)
    // Convert to human form
    lattice_simulator.convert_normal_to_human(particles);
}
