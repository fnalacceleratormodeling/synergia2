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

//void populate_6d_gaussian_quasi(Array_2d<double> &particles,
//                     const Array_1d<double> &means, const Array_2d<double> &covariances,
//                     const int id_offset);
//
//void populate_transverse_gaussian(Array_2d<double> &particles,
//                                  const Array_1d<double> &means, const Array_2d<double> &covariances,
//                                  const int id_offset, const unsigned long int seed, bool init_generator);
//
//void populate_transverse_gaussian_quasi(Array_2d<double> &particles,
//                                  const Array_1d<double> &means, const Array_2d<double> &covariances,
//                                  const int id_offset);
//
//void populate_uniform_cylinder(Array_2d<double> &particles,
//                             const Array_1d<double> &means, const Array_2d<double> &covariances,
//                             const int id_offset, const unsigned long int seed, bool init_generator);


#endif /* POPULATE_H_ */
