// -*- C++ -*-
#ifndef HAVE_FFTW2_HELPER_H
#define HAVE_FFTW2_HELPER_H true

#include "fftw_helper.h"
#include <rfftw_mpi.h>

class Fftw2_helper_nompi: public Fftw_helper
{
private:
  rfftwnd_plan plan, inv_plan;
  int upper_limit,max_local_size;
public:
  Fftw2_helper_nompi(Real_scalar_field &rho);
  int lower();
  int upper();
  size_t local_size();
  Int3 padded_shape_real();
  Int3 padded_shape_complex();
  void transform(Real_scalar_field &in, Complex_scalar_field &out);
  void inv_transform(Complex_scalar_field &in, Real_scalar_field &out);
  ~Fftw2_helper_nompi();
};

class Fftw2_helper_mpi: public Fftw_helper
{
private:
  rfftwnd_mpi_plan plan, inv_plan;
  int lower_limit,upper_limit,max_local_size;
  fftw_real *workspace, *data;
public:
  Fftw2_helper_mpi(Real_scalar_field &rho);
  int lower();
  int upper();
  size_t local_size();
  Int3 padded_shape_real();
  Int3 padded_shape_complex();
  void transform(Real_scalar_field &in, Complex_scalar_field &out);
  void inv_transform(Complex_scalar_field &in, Real_scalar_field &out);
  ~Fftw2_helper_mpi();
};

#endif
