#include "fftw3_helper.h"

Fftw3_helper_mpi::Fftw3_helper_mpi(Real_scalar_field &rho):
  Fftw_helper(rho)
{
  fftw_mpi_init();
  
  int local_n;
  Int3 shape(rho.get_points().get_shape());
  shape.scale(2);
  local_size = fftw_mpi_local_size(3,shape,MPI_COMM_WORLD,&local_n,
				    &local_lower);
  local_upper = local_lower + local_n;
  std::cout << "local size is " << local_size << ", would have guessed "
	    << local_n*shape[1]*shape[2] << std::endl;
  std::cout << "shape 0: " << shape[0] << ", local_lower = "
	    << local_lower << ", local_upper = " << local_upper << std::endl;

  a = (double *)fftw_malloc(local_size*sizeof(double));
  ahat = (fftw_complex *)fftw_malloc(local_size*
				     sizeof(fftw_complex));
  const int flags = FFTW_ESTIMATE;
  plan = fftw_mpi_plan_dft_r2c(3,shape.c_array(),a,ahat,MPI_COMM_WORLD,
			   flags);
  inv_plan = fftw_mpi_plan_dft_c2r(3,shape.c_array(),ahat,a,MPI_COMM_WORLD,
			       flags);
}

int 
Fftw3_helper_mpi::lower()
{
  return local_lower;
}

int 
Fftw3_helper_mpi::upper()
{
  return local_upper;
}

void
Fftw3_helper_mpi::transform(Real_scalar_field &in, Complex_scalar_field &out)
{
  memcpy(reinterpret_cast<void*>(a),
	 reinterpret_cast<void*>(in.get_points().get_base_address()),
	 in.get_points().get_length()*sizeof(double));
  std::cout << "copied " << in.get_points().get_length() 
	    << " doubles to a space for " << local_size << std::endl;
  fftw_execute(plan);
  memcpy(reinterpret_cast<void*>(out.get_points().get_base_address()),
	 reinterpret_cast<void*>(ahat),
	 out.get_points().get_length()*sizeof(std::complex<double>));
}

void 
Fftw3_helper_mpi::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
  memcpy(reinterpret_cast<void*>(ahat),
	 reinterpret_cast<void*>(in.get_points().get_base_address()),
	 in.get_points().get_length()*sizeof(std::complex<double>));
  fftw_execute(inv_plan);
  memcpy(reinterpret_cast<void*>(out.get_points().get_base_address()),
	 reinterpret_cast<void*>(a),
	 out.get_points().get_length()*sizeof(double));
}

Fftw3_helper_mpi::~Fftw3_helper_mpi()
{
  std::cout << "jfa: deleting a\n";
  fftw_free((void *)a);
  std::cout << "and ahat\n";
  fftw_free((void *)ahat);
  std::cout << "jfa: done deleting a and ahat\n";
}

Fftw3_helper_nompi::Fftw3_helper_nompi(Real_scalar_field &rho):
  Fftw_helper(rho)
{
  Int3 shape(rho.get_points().get_shape());
  shape.scale(2);
  upper_limit = shape[0];
  double *a;
  fftw_complex *ahat;
  a = (double *)fftw_malloc(shape[0]*shape[1]*shape[2]*sizeof(double));
  ahat = (fftw_complex *)fftw_malloc(shape[0]*shape[1]*(shape[2]/2+1)*
				     sizeof(fftw_complex));
  const int flags = FFTW_ESTIMATE;
  plan = fftw_plan_dft_r2c(3,shape.c_array(),a,ahat,
			   flags);
  inv_plan = fftw_plan_dft_c2r(3,shape.c_array(),ahat,a,
			       flags);
  fftw_free((void *)a);
  fftw_free((void *)ahat);
}

int 
Fftw3_helper_nompi::lower()
{
  return 0;
}

int 
Fftw3_helper_nompi::upper()
{
  return upper_limit;
}

void
Fftw3_helper_nompi::transform(Real_scalar_field &in, Complex_scalar_field &out)
{
  fftw_execute_dft_r2c(plan,
		       in.get_points().get_base_address(),
		       reinterpret_cast<fftw_complex *>
		       (out.get_points().get_base_address()));
}

void 
Fftw3_helper_nompi::inv_transform(Complex_scalar_field &in, Real_scalar_field &out)
{
  fftw_execute_dft_c2r(inv_plan,
		       reinterpret_cast<fftw_complex *>
		       (in.get_points().get_base_address()),
		       out.get_points().get_base_address());
}

Fftw3_helper_nompi::~Fftw3_helper_nompi()
{
}

