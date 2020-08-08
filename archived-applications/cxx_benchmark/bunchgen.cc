#include <iostream>
#include <stdexcept>

#include "synergia/foundation/reference_particle.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/core_diagnostics.h"

#include "synergia/utils/synergia_omp.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    const int num_macro_particles = 10485760;
    const int seed = 4;
    const double num_real_particles = 1e13;
    const double mass = 1.0;
    const double energy = 2.0;
    const double emittance = 12.0e-6;
    const double beta_x = 12.0;
    const double alfa_x = 0.1;
    const double beta_y = 33.0;
    const double alfa_y = -0.2;
    const double stdz = 0.5;
    const double stddpop= 1.0e-4;

    const Reference_particle refpart(1, mass, energy);

    Commxx_sptr comm_sptr(new Commxx);
    Bunch_sptr bunch_sptr(
            new Bunch(refpart,
                    num_macro_particles, num_real_particles, comm_sptr));
    Random_distribution distribution(seed, *comm_sptr);

    MArray1d means(boost::extents[6]);
    for (int i=0; i<6; ++i) {
        means[i] = 0.0;
    }

    MArray2d covariances(boost::extents[6][6]);
    for (int i=0; i<6; ++i) {
        for (int j=0; j<6; ++j) {
            covariances[i][j] = 0.0;
        }
    }

    covariances[0][0] = emittance * beta_x;
    covariances[0][1] = (-alfa_x/beta_x) * emittance * beta_x;
    covariances[1][0] = covariances[0][1];
    covariances[1][1] = emittance * (1.0 + alfa_x*alfa_x)/beta_x;
    covariances[2][2] = emittance * beta_y;
    covariances[2][3] = (-alfa_y/beta_y) * emittance * beta_y;
    covariances[3][2] = covariances[2][3];
    covariances[3][3] = emittance * (1.0 + alfa_y*alfa_y)/beta_y;
    covariances[4][4] = stdz*stdz;
    covariances[5][5] = stddpop*stddpop;

    if (comm_sptr->get_rank() == 0) {
        std::cout << "covariance matrix:" << std::endl;
        std::cout.precision(16);
        for (int i=0; i<6; ++i) {
            if (i != 0) {
                std::cout << std::endl;
            }
            for (int j=0; j<6; ++j) {
                if (j != 0) {
                    std::cout << " ";
                }
                std::cout << covariances[i][j];
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
    double t1 = MPI_Wtime();
    populate_6d(distribution, *bunch_sptr, means, covariances);
    double gentime = MPI_Wtime() - t1;
    MArray1d bunch_means(Core_diagnostics::calculate_mean(*bunch_sptr));
    MArray1d bunch_stds(Core_diagnostics::calculate_std(*bunch_sptr, bunch_means));
    if (comm_sptr->get_rank() == 0) {
        std::cout << "Bunch generation time: " << gentime << std::endl;
        std::cout.precision(16);
        std::cout << "Bunch means:";
        for (int i=0; i<6; ++i) {
            std::cout << " " << bunch_means[i];
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "Bunch stds:" << std::endl;
        bool bad = false;
        for (int i=0; i<6; ++i) {
            std::cout << i << ": " << bunch_stds[i] << "\t<-->  " << std::sqrt(covariances[i][i]) <<
                         "\tDiff: " << (bunch_stds[i] - std::sqrt(covariances[i][i]))/bunch_stds[i];
            if (abs((bunch_stds[i] - std::sqrt(covariances[i][i]))/bunch_stds[i]) > 1.0e-13) {
                std::cout << "  XXXXX!!!!!";
                bad = true;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        if (bad) {
            std::cout << "!!!!! Error! Bad bunch generation statistics!!!!" << std::endl;
        }
        std::cout << std::endl;
    }
}

int
main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
