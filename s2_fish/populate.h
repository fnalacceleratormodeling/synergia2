#ifndef HAVE_POPULATE_H
#define HAVE_POPULATE_H true

#include "array_1d.h"
#include "array_2d.h"

void populate_6d_gaussian(Array_2d<double> &particles, 
    const Array_1d<double> &means, const Array_2d<double> &covariances,
    const int id_offset);

#endif
