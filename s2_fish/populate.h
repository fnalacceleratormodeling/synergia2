#ifndef HAVE_POPULATE_H
#define HAVE_POPULATE_H true

#include "array_nd/array_1d.h"
#include "array_nd/array_2d.h"

void populate_6d_gaussian(Array_2d<double> &particles,
                          const Array_1d<double> &means, const Array_2d<double> &covariances,
                          const int id_offset, const unsigned long int seed, bool init_generator);

void populate_6d_gaussian_quasi(Array_2d<double> &particles,
                     const Array_1d<double> &means, const Array_2d<double> &covariances,
                     const int id_offset);

void populate_transverse_gaussian(Array_2d<double> &particles,
                                  const Array_1d<double> &means, const Array_2d<double> &covariances,
                                  const int id_offset, const unsigned long int seed, bool init_generator);

#endif
