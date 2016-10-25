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

/// Populate a bunch with a truncated Gaussian distribution in all six dimensions.
/// @param dist the distribution generator
/// @param bunch the bunch
/// @param means an array of length six of the mean value of the distribution
///  in each phase space variable
/// @param covariances the six-by-six covariance matrix
/// @param limits an array of length six giving the cutoffs in units of the individual sigmas. A zero value means "do not truncate."
void
populate_6d_truncated(Distribution &dist, Bunch &bunch,
        Const_MArray1d_ref means, Const_MArray2d_ref covariances,
        Const_MArray1d_ref limits);

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

/// Populate a bunch with a Kapichinskij and Vladimirskij (KV) distribution
/// on the transverse plane, flat in cdt and Gaussian in dpop.
/// Reference for KV distribution is J. D. Lawson, the Physics of 
/// charged particle beam, 2nd edition, p 170. 
/// @param dist the distribution generator for the Gaussian dp/p distribution. 
/// @param bunch the bunch
/// @param epsilMax_x: The maximum emittance (not normalized) to fill horizontal twoxtwo D ellipse. 
/// @param epsilMax_y: The maximum emittance (not normalized) to fill vertical twoxtwo D ellipse. 
/// @param aplha_x: Alpha Twiss parameter at the point of injection.
/// @param beta_x: Beta Twiss parameter at the point of injection.
/// @param aplha_y: Alpha Twiss parameter at the point of injection.
/// @param beta_y: Beta Twiss parameter at the point of injection.
/// @param cdt the total range of the longitudinal coordinate will be
///  [-cdt/2,cdt/2] [m] (flatly distributed). 
/// @param stddpop standard deviation of the dp/p distribution
void
populate_transverse_KV_GaussLong(Distribution &dist, Bunch &bunch, double epsilMax_x,
        double alpha_x, double beta_x, double epsilMax_y, double alpha_y, double beta_y,
        double stddt, double stddpop);
	
void 
populate_two_particles(Bunch &bunch, 
         double p1x, double p1xp, double p1y, double p1yp, double p1cdt, double p1dpop, 
         double p2x, double p2xp, double p2y, double p2yp, double p2cdt, double p2dpop); 

#endif /* POPULATE_H_ */
