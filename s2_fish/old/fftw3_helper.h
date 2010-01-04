// -*- C++ -*-
#ifndef HAVE_FFTW3_HELPER_H
#define HAVE_FFTW3_HELPER_H true

#include "scalar_field.h"
#include "fftw_helper.h"
#include <fftw3-mpi.h>

class Fftw3_helper_mpi: public Fftw_helper
{
private:
  fftw_plan plan, inv_plan;
  ptrdiff_t local_size, local_lower, local_upper;
  double *a;
  fftw_complex *ahat;
public:
  Fftw3_helper_mpi(Real_scalar_field &rho);
  int lower();
  int upper();
  void transform(Real_scalar_field &in, Complex_scalar_field &out);
  void inv_transform(Complex_scalar_field &in, Real_scalar_field &out);
  ~Fftw3_helper_mpi();
};

class Fftw3_helper_nompi: public Fftw_helper
{
private:
  fftw_plan plan, inv_plan;
  int upper_limit;
public:
  Fftw3_helper_nompi(Real_scalar_field &rho);
  int lower();
  int upper();
  void transform(Real_scalar_field &in, Complex_scalar_field &out);
  void inv_transform(Complex_scalar_field &in, Real_scalar_field &out);
  ~Fftw3_helper_nompi();
};

#endif
