#include "solver.h"
#include "petscvec.h"
#include "petscmat.h"
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>

#include "math_constants.h"

#include "petscda.h"
#include "petscksp.h"
#include "petscdmmg.h"

#include "mytimer.h"

void 
check(PetscErrorCode ierr)
{
  if (ierr != 0) {
    std::stringstream message("");
    message << "Petsc error in solver ierr = " << ierr;
    throw std::runtime_error(message.str());
  }
}

void
print_vec(std::string prefix, Real_scalar_field sf, bool doubled, Vec vec)
{
  PetscErrorCode ierr;
  Int3 shape(sf.get_points().get_shape());
  if (doubled) {
    shape.scale(2);
  }
  PetscScalar ***sf_array;
  ierr = VecGetArray3d(vec,shape[0],shape[1],shape[2],0,0,0,
                       &sf_array); check(ierr);
  Int3 point;
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

void
write_vec(std::string filename, Real_scalar_field sf, bool doubled, Vec vec)
{
  PetscErrorCode ierr;
  Int3 shape(sf.get_points().get_shape());
  if (doubled) {
    shape.scale(2);
  }
  PetscScalar *sf_array;
  ierr = VecGetArray1d(vec,shape[0]*shape[1]*shape[2],
		       0,&sf_array); check(ierr);
  std::ofstream stream(filename.c_str());
  stream << shape[0] << " "
	 << shape[1] << " "
	 << shape[2] << std::endl;
  for(int i=0; i<shape[0]*shape[1]*shape[2]; ++i) {
    stream << sf_array[i];
    stream << " ";
  }
  stream << std::endl;
  stream.close();
  ierr = VecRestoreArray1d(vec,shape[0],0,&sf_array); check(ierr);
}

double
fft_tester(int nx, int ny, int nz)
{
  int array_length = nx*ny*nz;
  int dims[3] = {nx,ny,nz};
  int argc = 0;
  char **argv;
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char *)0,""); check(ierr);
  PetscScalar data[array_length];
  for(int i=0; i<array_length; ++i) {
    data[i] = std::complex<double>(1.000001*i + 0.1,2.000001*i + 0.5);
  }
  Vec v;
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,array_length,
			       data,&v); check(ierr);
  Vec vhat;
  ierr = VecDuplicate(v,&vhat); check(ierr);

  Vec vhathat;
  ierr = VecDuplicate(v,&vhathat); check(ierr);

  Mat A;
  ierr = MatCreateSeqFFTW(PETSC_COMM_SELF,
			  3,dims,&A); check(ierr);

  ierr = MatMult(A,v,vhat); check(ierr);

  PetscScalar *aa,*aa1,*aa2;
  ierr = VecGetArray(v,&aa); check(ierr);
  ierr = VecGetArray(vhat,&aa1); check(ierr);

//   for(int i=0; i<array_length; ++i) {
//     std::cout << "in[" << i << "] = " << aa[i]
// 	      << ", out[" << i << "] = " << aa1[i] << std::endl;
//   }

  double norm = 1.0/(array_length);
  ierr = VecScale(vhat,norm);
  
  ierr = MatMultTranspose(A,vhat,vhathat); check(ierr);

  PetscScalar ***a,***a2;
  ierr = VecGetArray3d(v,dims[0],dims[1],dims[2],
		       0,0,0,&a); check(ierr);
  ierr = VecGetArray3d(vhathat,dims[0],dims[1],dims[2],
		       0,0,0,&a2); check(ierr);
  double max_err(0.0);
  double ratio;
  for(int i=0; i<dims[0]; ++i) {
    for(int j=0; j<dims[1]; ++j) {
      for(int k=0; k<dims[2]; ++k) {
	if (a2[i][j][k] != 0.0) {
	  ratio = a[i][j][k].real()/a2[i][j][k].real();
	} else {
	  ratio = a[i][j][k].real() + 1.0;
	}
	double abs_diff = fabs(ratio-1.0);
	if (abs_diff > max_err) {
	  max_err = abs_diff;
	}
      }
    }
  }
  return max_err;
}

Vec
get_rho_hat2_petsc(Real_scalar_field &rho, Mat &FFT_matrix)
{
  // steps 1 and 2
  Int3 num_points2 = rho.get_points().get_shape();
  num_points2.scale(2);
  Double3 physical_size2 = rho.get_physical_size();
  physical_size2.scale(2.0);
  Complex_scalar_field rho2(num_points2.vector(), physical_size2.vector(),
			    rho.get_physical_offset());
  rho2.get_points().zero_all();
  Int3 num_points = rho.get_points().get_shape();
  Int3 index;
  for(index[0]=0; index[0]<num_points[0]; ++index[0]) {
    for(index[1]=0; index[1]<num_points[1]; ++index[1]) {
      for(index[2]=0; index[2]<num_points[2]; ++index[2]) {
	rho2.get_points().set(index,rho.get_points().get(index));
      }
    }
  }
  Vec rho2_petsc;
  PetscErrorCode ierr;
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,rho2.get_points().get_length(),
			       rho2.get_points().get_base_address(),
			       &rho2_petsc); check(ierr);
  Vec rho_hat2_petsc;
  ierr = VecDuplicate(rho2_petsc,&rho_hat2_petsc); check(ierr);
  ierr = MatMult(FFT_matrix,rho2_petsc,
		 rho_hat2_petsc); check(ierr);
  ierr = VecDestroy(rho2_petsc); check(ierr);
  return rho_hat2_petsc;
}

Vec
get_G_hat2_petsc(Real_scalar_field &rho, Mat &FFT_matrix, bool z_periodic)
{
  //step 3
  const double pi = 4.0*atan(1.0);
  Int3 num_points2 = rho.get_points().get_shape();
  num_points2.scale(2);
  Double3 physical_size = rho.get_physical_size();
  Double3 physical_size2 = rho.get_physical_size();
  physical_size2.scale(2.0);
  Complex_scalar_field G2(num_points2.vector(), physical_size2.vector(),
			  rho.get_physical_offset());
  G2.get_points().zero_all();
  Double3 h(rho.get_cell_size());
  Int3 index;
  // What is G(0,0,0), anyway? Rob and Ji seem to think it is G(0,0,1).
  // Hockney seems to think it is 1.
  // I don't think it is either, but I have not yet worked out what I 
  // consider to be the right answer (the one that preserves the integral
  // of G*rho).
  
  // Rob and Ji version
//   G2.get_points().set(Int3(0,0,0),G2.get_points().get(Int3(0,0,1)));

  // This would be the correct value if we were using cells that were spheres
  // instead of rectangular solids:
  // average value of inner sphere:
//   G2.get_points().set(Int3(0,0,0),(1.0/4.0*pi)*
// 		      (3.0/(2.0*((1+sqrt(3))/2.0)*
// 			    sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]))));

  // average value of outer sphere. This works unreasonably well.
  double G000 = (1.0/4.0*pi)*(3.0/(2.0*(sqrt(3))*
				   sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2])));
  // Calculate what I think the answer should be by doing a very
  // simple numerical integral. The resulting answer doesn't work
  // nearly as well as the outer sphere approximation above, so
  // something must be wrong.
//   double sum = 0.0;
//   const int nsteps = 32;
//   double dx = h[0]/nsteps;
//   double dy = h[1]/nsteps;
//   double dz = h[2]/nsteps;
//   x = 0.5*dx;
//   for(int i=0; i<nsteps; ++i) {
//     y = 0.5*dy;
//     for(int j=0; j<nsteps; ++j) {
//       z = 0.5*dz;
//       for(int k=0; k<nsteps; ++k) {
// 	sum += 1.0/sqrt(x*x+y*y+z*z);
// 	z += dz;
//       }
//       y += dy;
//     }
//     x += dx;
//   }
// //   double R = h[0];
// //   double vol = (1/8.0)*4.0*pi/3.0*R*R*R;
//   double R = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
//   double vol = h[0]*h[1]*h[2];
//   double mean_inv_r = sum*dx*dy*dz/vol;
//   G2.get_points().set(Int3(0,0,0),(1.0/4.0*pi)*mean_inv_r);

  double x,y,z,G;
  const int num_images = 4;
  for(index[0]=0; index[0]<num_points2[0]; ++index[0]) {
    if(index[0]>=num_points2[0]/2) {
      x = (num_points2[0]-index[0])*h[0];
    } else {
      x = index[0]*h[0];
    }
    for(index[1]=0; index[1]<num_points2[1]; ++index[1]) {
      if(index[1]>=num_points2[1]/2) {
	y = (num_points2[1]-index[1])*h[1];
      } else {
	y = index[1]*h[1];
      }
      for(index[2]=0; index[2]<num_points2[2]; ++index[2]) {
	if(index[2]>=num_points2[2]/2) {
	  z = (num_points2[2]-index[2])*h[2];
	} else {
	  z = index[2]*h[2];
	}
	if (!((x==0.0) && (y==0.0) && (z==0.0))) {
	  G = 1.0/(4.0*pi*sqrt(x*x + y*y + z*z));
	} else {
	  G = G000;
	}
	if (z_periodic) {
	  for(int image=-num_images; image<=num_images; ++image) {
	    if (image != 0) {
	      double z_image = z+image*physical_size[2];
	      if (!((x==0.0) && (y==0.0) && (fabs(z_image)<1.0e-14))) {
		G += 1.0/(4.0*pi*sqrt(x*x + y*y + z_image*z_image));
	      } else {
		G += G000;
	      }
	    }
	  }
	}
	G2.get_points().set(index,G);
      }
    }
  }


  Vec G2_petsc;
  PetscErrorCode ierr;
  ierr = VecCreateSeqWithArray(PETSC_COMM_WORLD,G2.get_points().get_length(),
			       G2.get_points().get_base_address(),&G2_petsc); check(ierr);
  //  write_vec("G2_petsc",rho,true,G2_petsc);
  Vec G_hat2_petsc;
  ierr = VecDuplicate(G2_petsc,&G_hat2_petsc); check(ierr);
  ierr = MatMult(FFT_matrix,G2_petsc,G_hat2_petsc); check(ierr);
  ierr = VecDestroy(G2_petsc); check(ierr);
  return G_hat2_petsc;
}

Vec
get_phi_hat2_petsc(Real_scalar_field &rho, Vec &rho_hat2_petsc,
		   Vec &G_hat2_petsc)
{
  // step 4
  PetscErrorCode ierr;
  Vec phi_hat2_petsc;
  ierr = VecCreateSeq(PETSC_COMM_SELF,rho.get_points().get_length()*8,
		      &phi_hat2_petsc); check(ierr);
  ierr = VecZeroEntries(phi_hat2_petsc); check(ierr);
    
  Int3 shape(rho.get_points().get_shape());
  shape.scale(2);
  PetscScalar ***rho_hat2_array, ***G_hat2_array;
  ierr = VecGetArray3d(rho_hat2_petsc,shape[0],shape[1],shape[2],0,0,0,
                       &rho_hat2_array); check(ierr);
  ierr = VecGetArray3d(G_hat2_petsc,shape[0],shape[1],shape[2],0,0,0,
                       &G_hat2_array); check(ierr);
  PetscScalar ***phi_hat2_array;
  ierr = VecGetArray3d(phi_hat2_petsc,shape[0],shape[1],shape[2],0,0,0,
                       &phi_hat2_array); check(ierr);
  Double3 h(rho.get_cell_size());
  for(int i=0; i<shape[0]; ++i) {
    for(int j=0; j<shape[1]; ++j) {
      for(int k=0; k<shape[2]; ++k) {
	phi_hat2_array[i][j][k] = 	
	  rho_hat2_array[i][j][k]*
	  G_hat2_array[i][j][k]*h[0]*h[1]*h[2];
      }
    }
  }
  ierr = VecRestoreArray3d(rho_hat2_petsc,shape[0],shape[1],shape[2],0,0,0,
			   &rho_hat2_array); check(ierr);
  ierr = VecRestoreArray3d(G_hat2_petsc,shape[0],shape[1],shape[2],0,0,0,
			   &G_hat2_array); check(ierr);
  ierr = VecRestoreArray3d(phi_hat2_petsc,shape[0],shape[1],shape[2],0,0,0,
			   &phi_hat2_array); check(ierr);
  return phi_hat2_petsc;	
}

Vec
get_phi2_petsc(Real_scalar_field &rho, Vec &phi_hat2_petsc, Mat &FFT_matrix)
{
  // step 5
  PetscErrorCode ierr;
  Vec phi2_petsc;
  ierr = VecDuplicate(phi_hat2_petsc,&phi2_petsc); check(ierr);
  ierr = MatMultTranspose(FFT_matrix,phi_hat2_petsc,phi2_petsc); check(ierr);
  PetscInt size;
  ierr = VecGetSize(phi2_petsc,&size); check(ierr);
  double norm = 1.0/size;
  ierr = VecScale(phi2_petsc,norm);
  return phi2_petsc;
}

Real_scalar_field
get_phi(Real_scalar_field &rho, Vec &phi2_petsc)
{
  // step 6
  Real_scalar_field phi(rho.get_points().get_shape(),rho.get_physical_size(),
			rho.get_physical_offset());
  PetscErrorCode ierr;
  Int3 shape(phi.get_points().get_shape());
  Int3 shape2(shape);
  shape2.scale(2);
  PetscScalar ***phi2_array;
  ierr = VecGetArray3d(phi2_petsc,shape2[0],shape2[1],shape2[2],0,0,0,
                       &phi2_array); check(ierr);
  Int3 point;
  for(int i=0; i<shape[0]; ++i) {
    point[0] = i;
    for(int j=0; j<shape[1]; ++j) {
      point[1] = j;
      for(int k=0; k<shape[2]; ++k) {
	point[2] = k;
	phi.get_points().set(point,(phi2_array[i][j][k]).real());
      }
    }
  }
  ierr = VecRestoreArray3d(phi2_petsc,shape2[0],shape2[1],shape2[2],0,0,0,
			   &phi2_array); check(ierr);
  return phi;	
}


Real_scalar_field
solver_fft_open(Real_scalar_field &rho, bool z_periodic)
{
  // The plan: Solve del^2 phi = rho by:
  //  1) convert rho to rho2, where 2 suffix indicates
  //     that we are using Hockney's doubled grid. rho2 = zero when  
  //     index outside of boundaries of rho (and so on for other vars).
  //  2) get FFT(rho2) = rho_hat2
  //  3) get (Green function) G_hat2 (intrinsically complex)
  //  4) calculate phi_hat2 = rho_hat2 * G_hat2
  //  5) calculate phi2 = invFFT(phi_hat2)
  //  6) extract phi from phi2

  double t0,t1,t00 = time();
  int argc = 0;
  char **argv;
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char *)0,""); check(ierr);
  Mat FFT_matrix;
  Int3 num_points2 = rho.get_points().get_shape();
  num_points2.scale(2);
  ierr = MatCreateSeqFFTW(PETSC_COMM_SELF,
			  3,num_points2.c_array(),
			  &FFT_matrix); check(ierr);
  Vec rho_hat2_petsc = get_rho_hat2_petsc(rho,FFT_matrix);
  Vec G_hat2_petsc = get_G_hat2_petsc(rho,FFT_matrix,z_periodic);
  Vec phi_hat2_petsc = get_phi_hat2_petsc(rho, rho_hat2_petsc,
					  G_hat2_petsc);
  ierr = VecDestroy(rho_hat2_petsc); check(ierr);
  ierr = VecDestroy(G_hat2_petsc); check(ierr);
  Vec phi2_petsc = get_phi2_petsc(rho,phi_hat2_petsc,FFT_matrix);
  ierr = VecDestroy(phi_hat2_petsc); check(ierr);
  Real_scalar_field phi = get_phi(rho,phi2_petsc);
  ierr = VecDestroy(phi2_petsc); check(ierr);

  return phi;
}

Real_scalar_field
calculate_E_n_higher_order(Real_scalar_field &phi, int n)
{
  if((n<0) || (n>2)) {
    std::stringstream message("");
    message << "calculate_E_n: invalid argument n=" << n
	    <<". Argument be in range 0<=n<=2";
    throw std::invalid_argument(message.str());
  }
  Real_scalar_field E(phi.get_points().get_shape(),phi.get_physical_size(),
		      phi.get_physical_offset());
  Int3 shape(phi.get_points().get_shape());
  Double3 hi(phi.get_cell_size());
  double h(hi[n]),delta;
  Int3 point;
  Int3 offset_plus(0,0,0),offset_minus(0,0,0);
  offset_plus[n] = 1;
  offset_minus[n] = -1;
  double deriv;
  for(int i=0; i<shape[0]; ++i) {
    point[0] = i;
    for(int j=0; j<shape[1]; ++j) {
      point[1] = j;
      for(int k=0; k<shape[2]; ++k) {
	point[2] = k;
	if ((point[n] > 1) && (point[n] < shape[n] - 2)) {
	  Int3 p0(point),p1(point),p2(point),p3(point);
	  p0.add(offset_minus);
	  p0.add(offset_minus);
	  p1.add(offset_minus);
	  p2.add(offset_plus);
	  p3.add(offset_plus);
	  p3.add(offset_plus);
	  deriv = 1.0/(12.0*h) * 
	    ( - phi.get_points().get(p3) + 8.0*phi.get_points().get(p2)
	      - 8.0*phi.get_points().get(p1) + phi.get_points().get(p0) );
	} else {
	  Int3 p0(point),p1(point);
	  if (point[n] == 0) {
	    p1.add(offset_plus);
	    delta = h;
	  } else if(point[n] == shape[n] - 1) {
	    p0.add(offset_minus);
	    delta = h;
	  } else {
	    p0.add(offset_minus);
	    p1.add(offset_plus);
	    delta = 2.0*h;
	  }
	  deriv = (phi.get_points().get(p1) - phi.get_points().get(p0))/delta;
	}
	E.get_points().set(point,deriv);
      }
    }
  }
  return E;
}

Real_scalar_field
calculate_E_n(Real_scalar_field &phi, int n)
{
  if((n<0) || (n>2)) {
    std::stringstream message("");
    message << "calculate_E_n: invalid argument n=" << n
	    <<". Argument be in range 0<=n<=2";
    throw std::invalid_argument(message.str());
  }
  Real_scalar_field E(phi.get_points().get_shape(),phi.get_physical_size(),
		      phi.get_physical_offset());
  Int3 shape(phi.get_points().get_shape());
  Double3 hi(phi.get_cell_size());
  double h(hi[n]),delta;
  Int3 point;
  Int3 offset_plus(0,0,0),offset_minus(0,0,0);
  offset_plus[n] = 1;
  offset_minus[n] = -1;
  double deriv;
  for(int i=0; i<shape[0]; ++i) {
    point[0] = i;
    for(int j=0; j<shape[1]; ++j) {
      point[1] = j;
      for(int k=0; k<shape[2]; ++k) {
	point[2] = k;
	Int3 p0(point),p1(point);
	if (point[n] == 0) {
	  p1.add(offset_plus);
	  delta = h;
	} else if(point[n] == shape[n] - 1) {
	  p0.add(offset_minus);
	  delta = h;
	} else {
	  p0.add(offset_minus);
	  p1.add(offset_plus);
	  delta = 2.0*h;
	}
	deriv = (phi.get_points().get(p1) - phi.get_points().get(p0))/delta;
	E.get_points().set(point,deriv);
      }
    }
  }
  return E;
}

void
apply_E_n_kick(Real_scalar_field &E, int n_axis, double tau,
	       Macro_bunch_store &mbs)
{
  if((n_axis<0) || (n_axis>2)) {
    std::stringstream message("");
    message << "apply_E_n_kick: invalid argument n_axis=" << n_axis
	    <<". Argument be in range 0<=n_axis<=2";
    throw std::invalid_argument(message.str());
  }
  // jfa: I am taking this calculation of "factor" more-or-less
  // directly from Impact.  I should really redo it in terms that make
  // sense to me
  double gamma = -1* mbs.ref_particle(5);
  double beta = sqrt(gamma*gamma - 1.0)/gamma;
  const  double c = 2.99792458e8;
  const  double pi = 4.0*atan(1.0);

  // GACK! THIS IS HARDCODED FOR PROTONS! FIXME FIXME FIXME FIXME
  // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
  // FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME
  double mass = 0.93827231e9; // should I come from macro_bunch_store?
  double eps0 = 1.0/(4*pi*c*c*1.0e-7); // using c^2 = 1/(eps0 * mu0)
  double Brho = gamma*beta*mass/c;
  double perveance0 = mbs.total_current/(2*pi*eps0*Brho*gamma*gamma*\
					 beta*c*beta*c);
  double xk = mbs.units(0);

  double factor = pi*perveance0*gamma*beta*beta/xk * 4.0*pi;
  factor *= 1.0/mbs.total_num;
  if (n_axis == 2) {
    factor *= beta*gamma*gamma;
  } else {
    // (We think) this is for the Lorentz transformation of the transverse
    // E field.
    factor *= gamma;
  }
  int index = 2*n_axis+1; // for axis n_axis = (0,1,2) corresponding to x,y,z,
                          // in particle store indexing, px,py,pz = (1,3,5)
  double kick;
  for(int n=0; n<mbs.local_num; ++n) {
      kick = tau * factor * E.get_val(Double3(mbs.local_particles(0,n),
					      mbs.local_particles(2,n),
					      mbs.local_particles(4,n)));
    mbs.local_particles(index,n) -= kick;
  }
}

// BC types are not consistent, but match our needs:
// mixed is for x,y closed and z periodic
// pzneumann is for neumann (open) with the periodic z 
typedef enum {DIRICHLET, NEUMANN, MIXED, PZNEUMANN} BCType;

typedef struct {
  BCType        bcType;
  PetscInt gs;
  Real_scalar_field *rho;
} UserContext;

Real_scalar_field
get_real_scalar_field_from_vec(Real_scalar_field rho, Vec vec)
{
  Real_scalar_field retval(rho.get_points().get_shape(),
			   rho.get_physical_size(),
			   rho.get_physical_offset());
  PetscErrorCode ierr;
  Int3 shape(retval.get_points().get_shape());
  PetscScalar ***retval_array;
  ierr = VecGetArray3d(vec,shape[0],shape[1],shape[2],0,0,0,
                       &retval_array); check(ierr);
  Int3 point;
  for(int i=0; i<shape[0]; ++i) {
    point[0] = i;
    for(int j=0; j<shape[1]; ++j) {
      point[1] = j;
      for(int k=0; k<shape[2]; ++k) {
	point[2] = k;
	retval.get_points().set(point,(retval_array[i][j][k]).real());
      }
    }
  }
  ierr = VecRestoreArray3d(vec,shape[0],shape[1],shape[2],0,0,0,
			   &retval_array); check(ierr);
  return retval;	
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(DMMG dmmg,Vec b)
{
  PetscErrorCode ierr;
  PetscInt       mx,my,mz;
  PetscInt       xs,ys,zs;
  PetscInt       xm,ym,zm;
  PetscScalar    h;
  PetscInt size;
  PetscScalar    ***array;
  DA             da = (DA)dmmg->dm;
  
  PetscFunctionBegin;

  UserContext    *user = (UserContext *) dmmg->user;

  VecGetSize(b,&size);

  int rank;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = DAGetInfo((DA)dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DAVecGetArray(da, b, &array);CHKERRQ(ierr);
  Int3 point;
  for (int k=zs; k<zs+zm; k++){
    point[2] = k;
    for (int j=ys; j<ys+ym; j++){
      point[1] = j;
      for(int i=xs; i<xs+xm; i++){
	point[0] = i;
	array[i][j][k] = user->rho->get_points().get(point);
      }
    }
  }
  ierr = DAVecRestoreArray(da, b, &array);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  VecGetSize(b,&size);
  if (user->bcType != DIRICHLET) {
    MatNullSpace nullspace;
    ierr = KSPGetNullSpace(dmmg->ksp,&nullspace);CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullspace,b,PETSC_NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeJacobian"
PetscErrorCode ComputeJacobian(DMMG dmmg, Mat J, Mat jac)
{
  DA             da = (DA)dmmg->dm;
  UserContext    *user = (UserContext *) dmmg->user;
  PetscErrorCode ierr;
  PetscInt       num,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscInt       ileft,iright,jleft,jright,kleft,kright;
  PetscScalar    v[7],Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;
  MatStencil     row,col[7];

  PetscFunctionBegin;
  if (user->bcType != NEUMANN) {
    printf(" Boundary condition must be neumann\n");
    exit(1);
  }
  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0); check(ierr);

  Hx = user->rho->get_physical_size()[0]/ (PetscReal)(mx-1);
  Hy = user->rho->get_physical_size()[1]/ (PetscReal)(my-1);
  Hz = user->rho->get_physical_size()[2]/ (PetscReal)(mz-1);

  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); check(ierr);
  
  for (k=zs; k<zs+zm; k++){
    for (j=ys; j<ys+ym; j++){
      for(i=xs; i<xs+xm; i++){
        row.i = i; row.j = j; row.k = k;
	num = 0;
	v[0] = 2.0*(HxHydHz + HxHzdHy + HyHzdHx);
	col[0].i = i;col[0].j = j; col[0].k=k;
	num++;
	ileft = i-1;
	/* Neumann boundary conditions in x */
	if (ileft == -1) {
	  ileft = 0;
	  v[0] -= HyHzdHx;
	} else {
	  v[num] = -HyHzdHx; col[num].i = ileft; col[num].j = j; col[num].k=k;
	  num++;
	}
	iright = i+1;
	if (iright == mx) {
	  iright = mx-1;
	  v[0] -= HyHzdHx;
	} else {
	  v[num] = -HyHzdHx; col[num].i = iright; col[num].j = j; col[num].k=k;
	  num++;
	}
	jleft = j-1;
	/* Neumann boundary conditions in y */
	if (jleft == -1) {
	  jleft = 0;
	  v[0] -= HxHzdHy;
	} else {
	  v[num] = -HxHzdHy; col[num].i = i; col[num].j = jleft; col[num].k=k;
	  num++;
	}
	jright = j+1;
	if (jright == my) {
	  jright = my-1;
	  v[0] -= HxHzdHy;
	} else {
	  v[num] = -HxHzdHy; col[num].i = i; col[num].j = jright; col[num].k=k;
	  num++;
	}
	/* Neumann boundary conditions in z */
	kleft = k-1;
	if (kleft == -1) {
	  kleft = 0;
	  v[0] -= HxHydHz;
	} else {
	  v[num] = -HxHydHz; col[num].i = i; col[num].j = j; col[num].k=kleft;
	  num++;
	}
	kright = k+1;
	if (kright == mz) {
	  kright = mz-1;
	  v[0] -= HxHydHz;
	} else {
	  v[num] = -HxHydHz; col[num].i = i; col[num].j = j; col[num].k=kright;
	  num++;
	}
	ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES); check(ierr);
      }
    }
  }
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY); check(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY); check(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "solver_fd_multigrid_open"
Real_scalar_field
solver_fd_multigrid_open(Real_scalar_field &rho)
{
  int argc = 3;
  char fakeprog[] = "./solver_fd_multigrid_open";
  char precondflag1[] = "-mg_coarse_pc_type";
  char precondflag2[] = "cholesky";
  char * argv1[] = {fakeprog,precondflag1,precondflag2};
  char ** argv = argv1;
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char *)0,""); check(ierr);

  DMMG           *dmmg;
  PetscScalar    mone = -1.0;
  PetscReal      norm;
  UserContext    user;
  Vec            soln;
  DA             da;

  PetscFunctionBegin;
  const int dmmg_nlevels = 6;
  ierr = DMMGCreate(PETSC_COMM_WORLD,dmmg_nlevels,PETSC_NULL,&dmmg);
  check(ierr);
   // user stuff
  for (int l = 0; l < DMMGGetLevels(dmmg); l++) {
    ierr = DMMGSetUser(dmmg,l,&user); check(ierr);
  }
  user.bcType = NEUMANN;
  user.gs = 3;
  user.rho = &rho;
  if (user.bcType == PZNEUMANN) {
    ierr = DACreate3d(PETSC_COMM_WORLD,DA_ZPERIODIC,DA_STENCIL_STAR,
		      //		    3,3,9,
		      user.gs,user.gs,4*(user.gs -1) + 1,
		      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&da);
  }
  else { 
    ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
		      user.gs,user.gs,user.gs,
		      PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&da);
  }
  check(ierr); 

  ierr = DMMGSetDM(dmmg,(DM)da); check(ierr);
  ierr = DADestroy(da); check(ierr);

  ierr = DMMGSetKSP(dmmg,ComputeRHS,ComputeJacobian); check(ierr);

  //For all BC types with a Neumann component, do the majic
  if (user.bcType != DIRICHLET) {
    ierr = DMMGSetNullSpace(dmmg,PETSC_TRUE,0,PETSC_NULL); check(ierr);
  }

  ierr = DMMGSolve(dmmg); check(ierr);
  ierr = VecAssemblyBegin(DMMGGetx(dmmg)); check(ierr);
  ierr = VecAssemblyEnd(DMMGGetx(dmmg)); check(ierr);

  soln = DMMGGetx(dmmg);

  ierr = MatMult(DMMGGetJ(dmmg),DMMGGetx(dmmg),DMMGGetr(dmmg)); check(ierr);
  ierr = VecAXPY(DMMGGetr(dmmg),mone,DMMGGetRHS(dmmg)); check(ierr);
  ierr = VecNorm(DMMGGetr(dmmg),NORM_2,&norm); check(ierr);

  PetscInt size;
  VecGetSize(soln,&size);
  Real_scalar_field phi = get_real_scalar_field_from_vec(rho,soln);
  ierr = DMMGDestroy(dmmg); check(ierr);
  return phi;
}  
