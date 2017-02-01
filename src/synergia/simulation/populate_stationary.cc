#include "populate_stationary.h"

#include <iostream>
#include "synergia/bunch/populate.h"
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

    lattice_simulator.convert_normal_to_xyz(particles);
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
    for (int part = 0; part < num_part; ++part) {
        for (int c=0; c<3; ++c) {

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
    lattice_simulator.convert_normal_to_xyz(particles);    
    
}



void
populate_6d_stationary_gaussian_adjust(Distribution &dist, Bunch &bunch,
                const std::vector<double> actions,
                Lattice_simulator& lattice_simulator,  Const_MArray1d_ref means, Const_MArray2d_ref covariances)
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
    for (int part = 0; part < num_part; ++part) {
        for (int c=0; c<3; ++c) {

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
    lattice_simulator.convert_normal_to_xyz(particles);    
    adjust_moments(bunch, means, covariances);   
    
}

void
populate_6d_stationary_gaussian_truncated(Distribution &dist, Bunch &bunch,
                const std::vector<double> actions,
                Lattice_simulator& lattice_simulator,  Const_MArray1d_ref limits )
{
  
  // limit =xmax/sigma_x , xmax=aperture maximum radius 
  
    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];

    // fill even columns with uniform [0,1) which will be converted
    // to exponential distribution later
    for (int c=0; c<6; c+=2) {
        dist.fill_uniform(particles[boost::indices[range()][c]], 0.0, 1.);
    }

    // fill odd columns with random phases 0 - 2*pi
    for (int c=1; c<6; c+=2) {
        dist.fill_uniform(particles[boost::indices[range()][c]], 0.0, 2.0*mconstants::pi);
    }

    // go through and set the components in complex cartesian normal form space
    // eq. 13 of Theory and Praxis of Map Analysis
    // a_k = i \sqrt{I_k} e^{-i \phi_k }
    for (int part = 0; part < num_part; ++part) {
        for (int c=0; c<3; ++c) {

            // the action is drawn from an exponential distribution with mean
            // actions[c] according to -mean*log(1-r) where r is random number
            // [0,1).  We need the square root of the action for the normal form
            // coordinate according to Eq. 13.
            double action_cutoff=1.-exp(-0.5*limits[c]*limits[c]); 
            double square_root_action = sqrt(-actions[c]*log(1.0-particles[part][2*c]*action_cutoff));
            double phase = particles[part][2*c+1];
            particles[part][2*c] = square_root_action * sin(phase);
            particles[part][2*c+1] = -square_root_action * cos(phase);
        }
    }
    lattice_simulator.convert_normal_to_xyz(particles);    
    
}


// fill a 2d array with random actions and angles
void
fill_6d_normal_form_coords(Distribution &dist, MArray2d_ref nf_particles,
                           const std::vector<double> actions)
{
    fill_6d_normal_form_coords(dist, nf_particles[boost::indices[range()][range()]],
                               actions);
}

void
fill_6d_normal_form_coords(Distribution &dist, MArray2d_view nf_particles,
                           const std::vector<double> actions)
{
    int num_part = nf_particles.shape()[0];

    for (int c=0; c<6; c+=2) {
        // fill even columns with uniform [0,1) which will be converted
        // to exponential distribution later
        dist.fill_uniform(nf_particles[ boost::indices[range()][c] ], 0.0, 1.0);

        // fill odd columns with random phases 0 - 2*pi
        dist.fill_uniform(nf_particles[ boost::indices[range()][c+1] ], 0.0, 2.0*mconstants::pi);
    }

    // go through and set the components in complex cartesian normal form space
    // eq. 13 of Theory and Praxis of Map Analysis
    // a_k = i \sqrt{I_k} e^{-i \phi_k }
    for (int part = 0; part < num_part; ++part) {
        for (int c=0; c<3; ++c) {
            // the action is drawn from an exponential distribution with mean
            // actions[c] according to -mean*log(1-r) where r is random number
            // [0,1).  We need the square root of the action for the normal form
            // coordinate according to Eq. 13.
            double square_root_action = sqrt(-actions[c]*log(1.0-nf_particles[part][2*c]));
            double phase = nf_particles[part][2*c+1];
            nf_particles[part][2*c] = square_root_action * sin(phase);
            nf_particles[part][2*c+1] = -square_root_action * cos(phase);
        }
    }    
}

void
populate_6d_stationary_clipped_longitudinal_gaussian(Distribution &dist,
                                                     Bunch &bunch,
                                                     const std::vector<double> actions,
                                                     double cdt_min, double cdt_max,
                                                     Lattice_simulator& lattice_simulator)
{
    const int max_tries = 100;  // throw if I don't converge in this many tries

    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];

    for (int part=0; part<num_part; ++part) {

        MArray2d nf_particle(boost::extents[1][7]);
        MArray2d test_particle(boost::extents[1][7]);

        int curr_try = 0;
        while (curr_try < max_tries) {

            // generate the normal form coordinates.  I'm going to reuse
            // the same coordinates with 4 different phases to get a set
            // of particles that should meet the longitudinal cut
            fill_6d_normal_form_coords(dist, nf_particle, actions);

            bool good_particle = true;

            for (int phase=0; phase<4; ++phase) {
                test_particle[boost::indices[0][range(0,6)]] = nf_particle[boost::indices[0][range(0,6)]];

                // convert to xyz form and then check
                lattice_simulator.convert_normal_to_xyz(test_particle);

                if ((test_particle[0][Bunch::cdt] < cdt_min) ||
                    (test_particle[0][Bunch::cdt] > cdt_max)) {

                    good_particle = false;
                    break; // bad particle, quit now and try again
                }

                // good so far.  rotate phase and check again
                double a2r = nf_particle[0][4];
                double a2i = nf_particle[0][5];
                // imaginary -> real, -real -> imaginary
                nf_particle[0][4] = a2i;
                nf_particle[0][5] = -a2r;
            } // retest

            if (good_particle) {
                break; // all right!  I have a good one
            }

            ++curr_try;  // otherwise I'm going around to try a different normal form particle

        }

        if (curr_try == max_tries) {
            throw std::runtime_error("6d_populate_stationary_clipped_gaussian couldn't produce an acceptable particle");
        }
        // otherwise, this is a good one
        // this particle has been through three pi/2 rotations, but all the
        // phases are random anyway.
        particles[boost::indices[part][range(0,6)]] = test_particle[boost::indices[0][range(0,6)]];
    }
}


void
populate_6d_stationary_truncated_longitudinal_gaussian(Distribution &dist,
                                                       Bunch &bunch,
                                                       const std::vector<double> actions,
                                                       double n_sigma,
                                                       Lattice_simulator& lattice_simulator)
{
    MArray2d_ref particles(bunch.get_local_particles());
    int num_part = particles.shape()[0];

    for (int c=0; c<6; c+=2) {
        // fill even columns with uniform [0,1) which will be converted
        // to exponential distribution later
        dist.fill_uniform(particles[ boost::indices[range()][c] ], 0.0, 1.0);

        // fill odd columns with random phases 0 - 2*pi
        dist.fill_uniform(particles[ boost::indices[range()][c+1] ], 0.0, 2.0*mconstants::pi);
    }

    // go through and set the components in complex cartesian normal form space
    // eq. 13 of Theory and Praxis of Map Analysis
    // a_k = i \sqrt{I_k} e^{-i \phi_k }
    for (int part = 0; part < num_part;) {
        for (int c=0; c<3; ++c) {
            // the action is drawn from an exponential distribution with mean
            // actions[c] according to -mean*log(1-r) where r is random number
            // [0,1).  We need the square root of the action for the normal form
            // coordinate according to Eq. 13.
            double square_root_action = sqrt(-actions[c]*log(1.0-particles[part][2*c]));
            double phase = particles[part][2*c+1];
            particles[part][2*c] = square_root_action * sin(phase);
            particles[part][2*c+1] = -square_root_action * cos(phase);
        }
        // check whether we meet satisfy the longitudinal cut.  Check
        // absolute value of both real and imaginary to cover all
        // phases
        std::cout << "checking longitudinal cut particle " << part;
        if ((std::abs(particles[part][4]) >= n_sigma*sqrt(actions[2])) ||
            (std::abs(particles[part][5]) >= n_sigma*sqrt(actions[2]))) {
            // nope, have to regenerate this one
            std::cout << "  FAILED" << std::endl;
            dist.fill_uniform(particles[boost::indices[part][range(0,6,2)] ], 0.0, 1.0);
            dist.fill_uniform(particles[boost::indices[part][range(1,6,2)] ], 0.0, 2.0*mconstants::pi);
            // part not incremented
        } else {
            // yes, continue to next particle
            std::cout << " PASSED!" << std::endl;
            ++part;
        }
    }

    lattice_simulator.convert_normal_to_xyz(particles);
}
