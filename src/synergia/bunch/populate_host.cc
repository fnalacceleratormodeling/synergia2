
//#include "populate_host.h"

#include "Eigen/Eigen"
#include "Eigen/Cholesky"

using namespace Eigen;

void
adjust_moments_eigen( 
        double const * means,
        double const * covariances,
        double const * bunch_mean,
        double const * bunch_mom2,
        int num_particles,
        int num_particles_slots,
        double * particles )
{
    Matrix<double, 6, 6, Eigen::ColMajor> C(covariances);
    Matrix<double, 6, 6, Eigen::ColMajor> G(C.llt().matrixL());
    Matrix<double, 6, 6, Eigen::ColMajor> X(bunch_mom2);
    Matrix<double, 6, 6, Eigen::ColMajor> H(X.llt().matrixL());
    Matrix<double, 6, 6, Eigen::ColMajor> A(G * H.inverse());

    // jfa: dummy exists only to work around a bad interaction betwen
    //      Eigen3 and g++ 4.1.2
    std::stringstream dummy;
    dummy << C;

    Eigen::Map<Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
            rho7(particles, num_particles_slots, 7);

    Matrix<double, 1, 6 > rhobar6(bunch_mean);

    for (int part = 0; part < num_particles; ++part) {
        rho7.block<1, 6>(part, 0) -= rhobar6;
    }

    rho7.block(0, 0, num_particles, 6) *= A.transpose();

    Matrix<double, 1, 6> means6(means);

    for (int part = 0; part < num_particles; ++part) {
        rho7.block<1, 6>(part, 0) += means6;
    }
}

