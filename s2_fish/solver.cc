#include "solver.h"
#include "petscvec.h"
#include "petscmat.h"
#include <iomanip>
#include <cmath>
#include <string>

void 
check(PetscErrorCode ierr)
{
  if (ierr != 0) {
    std::stringstream message("");
    message << "Petsc error in solver ierr = " << ierr;
    throw std::runtime_error(message.str());
  }
}

double zero_cutoff = 1.0e-10;
// prints a 3d petsc_vector whose dimensions are that of scalar_field
void
print(Real_scalar_field scalar_field,Vec petsc_vector)
{  
  PetscScalar ***a;
  int3 dims(scalar_field.get_num_points());
  PetscErrorCode ierr = VecGetArray3d(petsc_vector,
				      dims[0],dims[1],dims[2],
				      0,0,0,&a); check(ierr);
  for(int k=0; k < dims[2]; ++k) {
    std::cout << "Re a[:][:][" << k << "]\n";
    for (int j=0; j<dims[1]; ++j) {
      for (int i=0; i < dims[0]; ++i) {
	std::cout << std::setw(9);
	PetscScalar tmp = a[i][j][k];
	if (fabs(tmp.real()) < zero_cutoff) {
	  tmp.real() = 0.0;
	}
	if (fabs(tmp.imag()) < zero_cutoff) {
	  tmp.imag() = 0.0;
	}
	std::cout << tmp.real();
      }
      std::cout << std::endl;
    }
  }

  for(int k=0; k < dims[2]; ++k) {
    std::cout << "Im a[:][:][" << k << "]\n";
    for (int j=0; j<dims[1]; ++j) {
      for (int i=0; i < dims[0]; ++i) {
	std::cout << std::setw(9);
	PetscScalar tmp = a[i][j][k];
	if (fabs(tmp.real()) < zero_cutoff) {
	  tmp.real() = 0.0;
	}
	if (fabs(tmp.imag()) < zero_cutoff) {
	  tmp.imag() = 0.0;
	}
	std::cout << tmp.imag();
      }
      std::cout << std::endl;
    }
  }
}

void
print_ratio(Real_scalar_field scalar_field,Vec petsc_vector, Vec petsc_vector2)
{  
  PetscScalar ***a,***a2;
  int3 dims(scalar_field.get_num_points());
  PetscErrorCode ierr = VecGetArray3d(petsc_vector,
				      dims[0],dims[1],dims[2],
				      0,0,0,&a); check(ierr);
  ierr = VecGetArray3d(petsc_vector2,
		       dims[0],dims[1],dims[2],
		       0,0,0,&a2); check(ierr);
  for(int k=0; k < dims[2]; ++k) {
    std::cout << "a[:][:][" << k << "]\n";
    for (int j=0; j<dims[1]; ++j) {
      for (int i=0; i < dims[0]; ++i) {
	std::cout << std::setw(9);
	PetscScalar tmp = a[i][j][k];
	if (fabs(tmp.real()) < zero_cutoff) {
	  tmp.real() = 0.0;
	}
	if (fabs(tmp.imag()) < zero_cutoff) {
	  tmp.imag() = 0.0;
	}
	PetscScalar tmp2 = a2[i][j][k];
	if (fabs(tmp2.real()) < zero_cutoff) {
	  tmp2.real() = 0.0;
	}
	if (fabs(tmp2.imag()) < zero_cutoff) {
	  tmp2.imag() = 0.0;
	}
	double ratio;
	if(tmp.real() == 0.0) {
	  ratio = 1.0;
	} else {
	  ratio = tmp.real()/tmp2.real();
	}
	std::cout << ratio;
      }
      std::cout << std::endl;
    }
  }
}

Real_scalar_field
fft_tester(Real_scalar_field rho)
{
  Real_scalar_field phi(rho.get_num_points(),rho.get_physical_size(),
			rho.get_physical_offset());
  phi.zero_the_points();
  int argc = 0;
  char **argv;
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char *)0,""); check(ierr);

  PetscScalar complex_rho[rho.array_length()];
  for(int i=0; i<rho.array_length(); ++i) {
    complex_rho[i] = *(rho.array_base_address()+i);
  }
  Vec v;
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,rho.array_length(),
			       complex_rho,&v); check(ierr);

  Vec vhat;
  ierr = VecDuplicate(v,&vhat); check(ierr);

  Vec vhathat;
  ierr = VecDuplicate(v,&vhathat); check(ierr);

  Mat A;
  ierr = MatCreateSeqFFTW(PETSC_COMM_SELF,
			  3,rho.get_num_points().c_array(),&A); check(ierr);

  ierr = MatMult(A,v,vhat); check(ierr);
  int3 dims = rho.get_num_points();
  double norm = 1.0/(dims[0]*dims[1]*dims[2]);
  ierr = VecScale(vhat,norm);
  
  /*  std::cout << "\nv =\n";
      print(rho,v);

      std::cout << "\nvhat =\n";
      print(rho,vhat);
  */
  ierr = MatMultTranspose(A,vhat,vhathat); check(ierr);

  /*  std::cout << "\nv =\n";
      print(rho,v);

      std::cout << "\nvhathat =\n";
      print(rho,vhathat);

      std::cout << "\nvhathat/v =\n";
      print_ratio(rho,vhathat,v);
  */
  return phi;
}

Vec
get_complex_rho_hat_petsc(Real_scalar_field rho, Mat FFT_matrix)
{
  // steps 1 and 2
  //old way
  //  PetscScalar complex_rho[rho.array_length()];
  //  for(int i=0; i<rho.array_length(); ++i) {
  //    complex_rho[i] = *(rho.array_base_address()+i);
  //  }
  //  Vec complex_rho_petsc;
  //  PetscErrorCode ierr;
  //  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,rho.array_length(),
  //			       complex_rho,&complex_rho_petsc); check(ierr);
  Complex_scalar_field complex_rho;
  complex_rho.copy(&rho);
  std::cout << "complex_rho = " << std::endl;
  complex_rho.get_points().print();
  Vec complex_rho_petsc;
  PetscErrorCode ierr;
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,complex_rho.array_length(),
			       complex_rho.array_base_address(),
			       &complex_rho_petsc); check(ierr);
  Vec complex_rho_hat_petsc;
  ierr = VecDuplicate(complex_rho_petsc,&complex_rho_hat_petsc); check(ierr);
  ierr = MatMult(FFT_matrix,complex_rho_petsc,
		 complex_rho_hat_petsc); check(ierr);
  ierr = VecDestroy(complex_rho_petsc); check(ierr);
  return complex_rho_hat_petsc;
}

Vec
get_G_hat_petsc(Real_scalar_field rho, Mat FFT_matrix)
{
  //step 3
  std::vector<int> shape(rho.get_points().get_shape());
  Complex_scalar_field G(rho.get_num_points(),rho.get_physical_size(),
			 rho.get_physical_offset());
  double3 h(rho.get_cell_size());
  int3 point;
  for(int i=0; i<shape[0]; ++i) {
    double x = i*h[0];
    point[0] = i;
    for(int j=0; j<shape[1]; ++j) {
      double y = j*h[1];
      point[1] = j;
      for(int k=0; k<shape[2]; ++k) {
	double z = k*h[2];
	point[2] = k;
	if (!((i==0) && (j==0) && (k==0))) {
	  G.set_point(point,1.0/sqrt(x*x + y*y + z*z));
	} else {
	  G.set_point(point,1.0);
	}
      }
    }
  }
  Vec G_petsc;
  PetscErrorCode ierr;
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,G.array_length(),
			       G.array_base_address(),&G_petsc); check(ierr);
  Vec G_hat_petsc;
  ierr = VecDuplicate(G_petsc,&G_hat_petsc); check(ierr);
  ierr = MatMult(FFT_matrix,G_petsc,G_hat_petsc); check(ierr);
  ierr = VecDestroy(G_petsc); check(ierr);
  return G_hat_petsc;
}

Vec
get_phi_hat2_petsc(Real_scalar_field rho, Vec complex_rho_hat_petsc,
		   Vec G_hat_petsc)
{
  // step 4
  PetscErrorCode ierr;
  Vec phi_hat2_petsc;
  ierr = VecCreateSeq(PETSC_COMM_SELF,rho.array_length()*8,
		      &phi_hat2_petsc); check(ierr);
  ierr = VecZeroEntries(phi_hat2_petsc); check(ierr);
    
  int3 shape(rho.get_points().get_shape());
  PetscScalar ***complex_rho_hat_array, ***G_hat_array;
  ierr = VecGetArray3d(complex_rho_hat_petsc,shape[0],shape[1],shape[2],0,0,0,
                       &complex_rho_hat_array); check(ierr);
  ierr = VecGetArray3d(G_hat_petsc,shape[0],shape[1],shape[2],0,0,0,
                       &G_hat_array); check(ierr);
  int3 shape2;
  shape2[0] = shape[0]*2;
  shape2[1] = shape[1]*2;
  shape2[2] = shape[2]*2;
  PetscScalar ***phi_hat2_array;
  ierr = VecGetArray3d(phi_hat2_petsc,shape2[0],shape2[1],shape2[2],0,0,0,
                       &phi_hat2_array); check(ierr);
  for(int i=0; i<shape[0]; ++i) {
    for(int j=0; j<shape[1]; ++j) {
      for(int k=0; k<shape[2]; ++k) {
	phi_hat2_array[i][j][k] = 	
	  complex_rho_hat_array[i][j][k]*
	  G_hat_array[i][j][k];
      }
    }
  }
  ierr = VecRestoreArray3d(complex_rho_hat_petsc,shape[0],shape[1],shape[2],0,0,0,
			   &complex_rho_hat_array); check(ierr);
  ierr = VecRestoreArray3d(G_hat_petsc,shape[0],shape[1],shape[2],0,0,0,
			   &G_hat_array); check(ierr);
  ierr = VecRestoreArray3d(phi_hat2_petsc,shape2[0],shape2[1],shape2[2],0,0,0,
			   &phi_hat2_array); check(ierr);
  return phi_hat2_petsc;	
}

Vec
get_phi2_petsc(Real_scalar_field rho, Vec phi_hat2_petsc)
{
  // step 5
  PetscErrorCode ierr;
  Mat FFT2_matrix;
  PetscInt *shape = rho.get_num_points().c_array();
  shape[0] *= 2;
  shape[1] *= 2;
  shape[2] *= 2;
  ierr = MatCreateSeqFFTW(PETSC_COMM_SELF,
			  3,shape,
			  &FFT2_matrix); check(ierr);
  Vec phi2_petsc;
  ierr = VecDuplicate(phi_hat2_petsc,&phi2_petsc); check(ierr);
  ierr = MatMultTranspose(FFT2_matrix,phi_hat2_petsc,phi2_petsc); check(ierr);
  PetscInt size;
  ierr = VecGetSize(phi2_petsc,&size); check(ierr);
  double norm = 1.0/size;
  ierr = VecScale(phi2_petsc,norm);
  return phi2_petsc;
}

Real_scalar_field
get_phi(Real_scalar_field rho, Vec phi2_petsc)
{
  // step 6
  Real_scalar_field phi(rho.get_num_points(),rho.get_physical_size(),
			rho.get_physical_offset());
  PetscErrorCode ierr;
  int3 shape(phi.get_points().get_shape());
  int3 shape2;
  shape2[0] = shape[0]*2;
  shape2[1] = shape[1]*2;
  shape2[2] = shape[2]*2;
  PetscScalar ***phi2_array;
  ierr = VecGetArray3d(phi2_petsc,shape2[0],shape2[1],shape2[2],0,0,0,
                       &phi2_array); check(ierr);
  int3 point;
  for(int i=0; i<shape[0]; ++i) {
    point[0] = i;
    for(int j=0; j<shape[1]; ++j) {
      point[1] = j;
      for(int k=0; k<shape[2]; ++k) {
	point[2] = k;
	phi.set_point(point,(phi2_array[i][j][k]).real());
      }
    }
  }
  ierr = VecRestoreArray3d(phi2_petsc,shape2[0],shape2[1],shape2[2],0,0,0,
			   &phi2_array); check(ierr);
  return phi;	
}

void
print_vec(std::string prefix, Real_scalar_field sf, bool doubled, Vec vec)
{
  PetscErrorCode ierr;
  int3 shape(sf.get_points().get_shape());
  if (doubled) {
    shape[0] = shape[0]*2;
    shape[1] = shape[1]*2;
    shape[2] = shape[2]*2;
  }
  PetscScalar ***sf_array;
  ierr = VecGetArray3d(vec,shape[0],shape[1],shape[2],0,0,0,
                       &sf_array); check(ierr);
  int3 point;
  for(int k=0; k<shape[0]; ++k) {
    std::cout << prefix << "(:,:," << k << ")\n";
    for(int j=0; j<shape[1]; ++j) {
      for(int i=0; i<shape[2]; ++i) {
	std::cout << sf_array[i][j][k];
      std::cout << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  ierr = VecRestoreArray3d(vec,shape[0],shape[1],shape[2],0,0,0,
			   &sf_array); check(ierr);
}

Real_scalar_field
solver(Real_scalar_field rho)
{
  // The plan: Solve del^2 phi = rho by:
  //  1) convert rho to complex_rho
  //  2) get FFT(complex_rho) = complex_rho_hat
  //  3) get (Green function) G_hat (intrinsically complex)
  //  4) calculate phi_hat2 = complex_rho_hat * G_hat, where 2 suffix indicates
  //     that we are using Hockney's doubled grid. phi_hat2 = zero when index outside
  //     of boundaries of complex_rho_hat (or G_hat).
  //  5) calculate phi2 = invFFT(phi_hat2)
  //  6) extract phi from phi2

  std::cout << "rho:" << std::endl;
  rho.get_points().print();

  int argc = 0;
  char **argv;
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char *)0,""); check(ierr);

  Mat FFT_matrix;
  ierr = MatCreateSeqFFTW(PETSC_COMM_SELF,
			  3,rho.get_num_points().c_array(),
			  &FFT_matrix); check(ierr);
  Vec complex_rho_hat_petsc = get_complex_rho_hat_petsc(rho,FFT_matrix);
  print_vec("complex_rho_hat_petsc",rho,false,complex_rho_hat_petsc);
  Vec G_hat_petsc = get_G_hat_petsc(rho,FFT_matrix);
  print_vec("G_hat_petsc",rho,false,G_hat_petsc);
  Vec phi_hat2_petsc = get_phi_hat2_petsc(rho, complex_rho_hat_petsc,
					  G_hat_petsc);
  ierr = VecDestroy(complex_rho_hat_petsc); check(ierr);
  ierr = VecDestroy(G_hat_petsc); check(ierr);
  Vec phi2_petsc = get_phi2_petsc(rho,phi_hat2_petsc);
  ierr = VecDestroy(phi_hat2_petsc); check(ierr);
  Real_scalar_field phi = get_phi(rho,phi2_petsc);
  ierr = VecDestroy(phi2_petsc); check(ierr);
  return phi;
}
