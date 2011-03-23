#ifndef POPULATE_H_
#define POPULATE_H_

#include "synergia/foundation/distribution.h"
#include "synergia/bunch/bunch.h"

/// Populate a bunch with a Gaussian distribution in all six dimensions.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param means an array of length six of the mean value of the distribution
///  in each phase space variable
/// @param covariances the six-by-six covariance matrix
void
populate_6d(Distribution &dist, Bunch &bunch, Const_MArray1d_ref means,
        Const_MArray2d_ref covariances);

/// Populate a bunch with a Gaussian distribution in all four
/// transverse dimensions. The time distribution is uniform, but the
/// dp/p distribution is also Gaussian.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param means an array of length six of the mean value of the distribution
///  in each phase space variable
/// @param covariances the six-by-six covariance matrix
/// @param cdt the total range of the longitudinal coordinate will be
///  [-cdt/2,cdt/2] [m]
void
populate_transverse_gaussian(Distribution &dist, Bunch &bunch,
        Const_MArray1d_ref means, Const_MArray2d_ref covariances, double cdt);

/// Populate a bunch with a uniform cylindrical spatial distribution. The
/// momentum-like variables are Gaussian distributed and uncorrelated.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param radius radius in of the cylinder in the x-y plane [m]
/// @param cdt the total range of the longitudinal coordinate will be
///  [-cdt/2,cdt/2] [m]
/// @param stdxp standard deviation of the xp distribution
/// @param stdyp standard deviation of the yp distribution
/// @param stddpop standard deviation of the dp/p distribution
void
populate_uniform_cylinder(Distribution &dist, Bunch &bunch, double radius,
        double cdt, double stdxp, double stdyp, double stddpop);

#endif /* POPULATE_H_ */
