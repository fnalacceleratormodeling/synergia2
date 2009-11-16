#include "populate.h"
#include "diagnostics.h"

#include "utils/eigen2/Eigen/Eigen"
#include "utils/eigen2/Eigen/Cholesky"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

#include "utils/multi_array_print.h"

void
adjust_moments(Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    Diagnostics_full2 diagnostics(bunch, 0.0);
    Matrix<double, 6, 6, Eigen::RowMajor > C(covariances.origin());
    Matrix<double, 6, 6, Eigen::RowMajor > X(diagnostics.get_mom2().origin());
    Matrix<double, 6, 6, Eigen::RowMajor > A = C.llt().matrixL() * X.llt().matrixL().inverse();

    // The suffix 7 is used in the following section to indicate that we
    // are using the 7-dimensional data structure with the real 6-dimensional
    // data plus the id (seventh dimension).
    Matrix<double, 7, 7 > A7;
    A7 << A, MatrixXd::Zero(6, 1), MatrixXd::Zero(1, 6), MatrixXd::Identity(1,
            1);
    Eigen::Map<Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >
            rho7(bunch.get_local_particles().origin(), bunch.get_local_num(), 7);
    Matrix<double,1,6> rhobar6(diagnostics.get_mean().origin());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        rho7.block<1,6>(part,0) -= rhobar6;
    }

    rho7 *= A7.transpose();

    Matrix<double,1,6> means6(means.origin());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        rho7.block<1,6>(part,0) += means6;
    }
}

void
populate_6d(Distribution &dist, Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    MArray2d_ref particles(bunch.get_local_particles());
    for (int i = 0; i < 6; ++i) {
        dist.fill_unit_gaussian(particles[boost::indices[range()][i]]);
    }
    adjust_moments(bunch, means, covariances);
}
