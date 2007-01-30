#include "solver.h"
#include "petscvec.h"
#include "petscmat.h"
#include <iomanip>
#include <cmath>

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
print(Scalar_Field scalar_field,Vec petsc_vector)
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
print_ratio(Scalar_Field scalar_field,Vec petsc_vector, Vec petsc_vector2)
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

Scalar_Field
solver(Scalar_Field rho)
{
  // solves I phi = 0
  Scalar_Field phi(rho.get_num_points(),rho.get_physical_size(),
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

  Vec vtilde;
  ierr = VecDuplicate(v,&vtilde); check(ierr);

  Vec vtildetilde;
  ierr = VecDuplicate(v,&vtildetilde); check(ierr);

  Mat A;
  ierr = MatCreateSeqFFTW(PETSC_COMM_SELF,
			  3,rho.get_num_points().c_array(),&A); check(ierr);

  ierr = MatMult(A,v,vtilde); check(ierr);

  std::cout << "\nv =\n";
  print(rho,v);

  std::cout << "\nvtilde =\n";
  print(rho,vtilde);

  ierr = MatMultTranspose(A,vtilde,vtildetilde); check(ierr);

  std::cout << "\nv =\n";
  print(rho,v);

  std::cout << "\nvtildetilde =\n";
  print(rho,vtildetilde);

  std::cout << "\nvtildetilde/v =\n";
  print_ratio(rho,vtildetilde,v);

  return phi;
}
