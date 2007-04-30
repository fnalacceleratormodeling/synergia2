// -*- C++ -*-
#ifndef HAVE_FFTW_HELPER_H
#define HAVE_FFTW_HELPER_H true

#include "scalar_field.h"

class Fftw_helper
{
private:
public:
    Int3 shape;
    Fftw_helper(Real_scalar_field &rho);
    virtual int lower();
    virtual int upper();
    virtual int guard_lower();
    virtual int guard_upper();
    virtual int offset();
    virtual size_t local_size();
    virtual Int3 padded_shape_real();
    virtual Int3 padded_shape_complex();
    virtual void transform(Real_scalar_field &in, Complex_scalar_field &out);
    virtual void inv_transform(Complex_scalar_field &in, Real_scalar_field &out);
    ~Fftw_helper();
};

#endif
