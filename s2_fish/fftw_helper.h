// -*- C++ -*-
#ifndef HAVE_FFTW_HELPER_H
#define HAVE_FFTW_HELPER_H true

#ifdef USE_FFTW2
#include <rfftw_mpi.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif //USE_FFTW2

#include "scalar_field.h"
#undef DL_IMPORT
#include <vector>


class Fftw_helper
{
private:
#ifdef USE_FFTW2
    rfftwnd_mpi_plan plan, inv_plan;
    fftw_real *workspace, *data;
#else
    fftw_plan plan, inv_plan;
    double *data;
    fftw_complex *workspace;
#endif //USE_FFTW2
    int lower_limit, upper_limit, max_local_size;
    int left_guard, right_guard;
    Int3 shape;
    void construct(int *shape, bool z_periodic);
    bool have_local_data;
public:
    Fftw_helper(Int3 shape, bool z_periodic);
    Fftw_helper(std::vector<int> shape, bool z_periodic);
    int lower();
    int upper();
    int guard_lower();
    int guard_upper();
    int offset();
    size_t local_size();
    Int3 padded_shape_real();
    Int3 padded_shape_complex();
    void transform(Real_scalar_field &in, Complex_scalar_field &out);
    void inv_transform(Complex_scalar_field &in, Real_scalar_field &out);
    ~Fftw_helper();
};

#endif
