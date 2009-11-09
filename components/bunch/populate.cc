#include "populate.h"
#include "diagnostics.h"

#include "utils/eigen2/Eigen/Core"
#include "utils/eigen2/Eigen/LU"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

void
adjust_moments(Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances)
{
    MArray2d scratch(boost::extents[bunch.get_local_num()][6]);
    Diagnostics_full2 diagnostics(bunch,0.0);
//    MArray1d actual_means(diagnostics.get_mean());
//    MArray2d actual_covs(diagnostics.get_mom2());
    Matrix<double, 6, 6 > C(diagnostics.get_mom2().origin());
//    Matrix<double, 6,6 > G(C.llt().compute());


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
