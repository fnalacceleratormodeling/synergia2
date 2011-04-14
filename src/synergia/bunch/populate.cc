#include "populate.h"
#include "diagnostics.h"

#include "synergia/utils/eigen2/Eigen/Eigen"
#include "synergia/utils/eigen2/Eigen/Cholesky"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

#include "synergia/utils/multi_array_print.h"

void
adjust_moments(Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    MArray1d bunch_mean(Diagnostics::calculate_mean(bunch));
    MArray2d bunch_mom2(Diagnostics::calculate_mom2(bunch, bunch_mean));
    Matrix<double, 6, 6, Eigen::RowMajor > C(covariances.origin());
    Matrix<double, 6, 6, Eigen::RowMajor > X(bunch_mom2.origin());
    Matrix<double, 6, 6, Eigen::RowMajor > A = C.llt().matrixL()
            * X.llt().matrixL().inverse();

    int num_particles = bunch.get_local_num();
    Eigen::Map<Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >
            rho7(bunch.get_local_particles().origin(), num_particles, 7);
    Matrix<double, 1, 6 > rhobar6(bunch_mean.origin());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        rho7.block<1, 6 > (part, 0) -= rhobar6;
    }

    rho7.block(0, 0, num_particles, 6) *= A.transpose();

    Matrix<double, 1, 6 > means6(means.origin());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        rho7.block<1, 6 > (part, 0) += means6;
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

void
populate_transverse_gaussian(Distribution &dist, Bunch &bunch,
        Const_MArray1d_ref means, Const_MArray2d_ref covariances, double cdt)
{
    MArray2d_ref particles(bunch.get_local_particles());
    for (int i = 0; i < 4; ++i) {
        dist.fill_unit_gaussian(particles[boost::indices[range()][i]]);
    }
    dist.fill_uniform(particles[boost::indices[range()][4]], 0.0, 1.0);
    dist.fill_unit_gaussian(particles[boost::indices[range()][5]]);

    MArray1d means_modified(means);
    means_modified[4] = 0.0;

    // Symmetry requires no correlations with the cdt coordinate. Make a copy
    // of the covariance matrix and manually set all correlations to zero.
    MArray2d covariances_modified(covariances);
    for (int k = 0; k < 6; ++k) {
        covariances_modified[k][4] = covariances_modified[4][k] = 0.0;
    }
    covariances_modified[4][4] = cdt * cdt / 12.0;
    adjust_moments(bunch, means_modified, covariances_modified);
}

void
populate_uniform_cylinder(Distribution &dist, Bunch &bunch, double radius,
        double cdt, double stdxp, double stdyp, double stddpop)
{
    MArray2d_ref particles(bunch.get_local_particles());
    dist.fill_unit_disk(particles[boost::indices[range()][Bunch::x]],
            particles[boost::indices[range()][Bunch::y]]);
    dist.fill_uniform(particles[boost::indices[range()][Bunch::cdt]], -cdt
            / 2.0, cdt / 2.0);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::xp]]);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::yp]]);
    dist.fill_unit_gaussian(particles[boost::indices[range()][Bunch::dpop]]);

    for (int part = 0; part < bunch.get_local_num(); ++part) {
        particles[part][Bunch::x] *= radius;
        particles[part][Bunch::y] *= radius;
        particles[part][Bunch::xp] *= stdxp;
        particles[part][Bunch::yp] *= stdyp;
        particles[part][Bunch::dpop] *= stddpop;
    }
}
