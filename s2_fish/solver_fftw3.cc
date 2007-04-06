#include "solver.h"
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>

#undef DL_IMPORT
#include <fftw3.h>

#include "mytimer.h"

namespace solver_fftw3 {
  Complex_scalar_field
  get_rho_hat2(Real_scalar_field &rho, fftw_plan &plan)
  {
    // steps 1 and 2
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field rho2(num_points2.vector(), physical_size2.vector(),
			   rho.get_physical_offset());
    timer("misc");
    rho2.get_points().zero_all();
    timer("calc zero all rho2");
    Int3 num_points = rho.get_points().get_shape();
    Int3 index;
    timer("misc");
    for(index[0]=0; index[0]<num_points[0]; ++index[0]) {
      for(index[1]=0; index[1]<num_points[1]; ++index[1]) {
	for(index[2]=0; index[2]<num_points[2]; ++index[2]) {
	  rho2.get_points().set(index,rho.get_points().get(index));
	}
      }
    }
    timer("calc rho2");
    Int3 complex_num_points(num_points2);
    complex_num_points[2] = complex_num_points[2]/2 + 1;
    Complex_scalar_field rho_hat2(complex_num_points.vector(),
				  physical_size2.vector(),
				  rho.get_physical_offset());
    timer("rho plan");
    fftw_execute_dft_r2c(plan,
			rho2.get_points().get_base_address(),
			reinterpret_cast<fftw_complex *>
			(rho_hat2.get_points().get_base_address()));
    timer("rho fft");
    return rho_hat2;
  }

  Complex_scalar_field
  get_G_hat2(Real_scalar_field &rho, bool z_periodic, fftw_plan plan)
  {
    //step 3
    const double pi = 4.0*atan(1.0);
    Int3 num_points = rho.get_points().get_shape();
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
    Double3 physical_size = rho.get_physical_size();
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field G2(num_points2.vector(), physical_size2.vector(),
			 rho.get_physical_offset());
    //jfa: not necessary!    G2.get_points().zero_all();
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
    timer("misc");
    double x,y,z,G;
    const int num_images = 4;
    int mix,miy,miz; // mirror index x, etc.
    for(index[0]=0; index[0]<=num_points[0]; ++index[0]) {
      x = index[0]*h[0];
      mix = num_points2[0] - index[0];
      if (mix == num_points[0]) {
	mix = num_points2[0];
      }
      for(index[1]=0; index[1]<=num_points[1]; ++index[1]) {
	y = index[1]*h[1];
	miy = num_points2[1] - index[1];
	if (miy == num_points[1]) {
	  miy = num_points2[1];
	}
	for(index[2]=0; index[2]<=num_points[2]; ++index[2]) {
	  z = index[2]*h[2];
	  miz = num_points2[2] - index[2];
	  if (miz == num_points[2]) {
	    miz = num_points2[2];
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
	  // seven mirror images
	  if (mix < num_points2[0]) {
	    G2.get_points().set(Int3(mix,index[1],index[2]),G);
	    if (miy < num_points2[1]) {
	      G2.get_points().set(Int3(mix,miy,index[2]),G);
	      if (miz < num_points2[2]) {
	      G2.get_points().set(Int3(mix,miy,miz),G);
	      }
	    }
	    if (miz < num_points2[2]) {
	      G2.get_points().set(Int3(mix,index[1],miz),G);
	    }
	  }
	  if (miy < num_points2[1]) {
	    G2.get_points().set(Int3(index[0],miy,index[2]),G);
	    if (miz < num_points2[2]) {
	      G2.get_points().set(Int3(index[0],miy,miz),G);
	    }
	  }
	  if (miz < num_points2[2]) {
	    G2.get_points().set(Int3(index[0],index[1],miz),G);
	  }
	}
      }
    }

    timer("calc G");
    Int3 complex_num_points(num_points2);
    complex_num_points[2] = complex_num_points[2]/2 + 1;
    Complex_scalar_field G_hat2(complex_num_points.vector(),
				physical_size2.vector(),
				rho.get_physical_offset());
    fftw_execute_dft_r2c(plan,
			 G2.get_points().get_base_address(),
			 reinterpret_cast<fftw_complex *>
			 (G_hat2.get_points().get_base_address()));
    timer("G fft");
    return G_hat2;
  }

  Complex_scalar_field
  get_phi_hat2(Real_scalar_field &rho, Complex_scalar_field &rho_hat2,
	       Complex_scalar_field &G_hat2)
  {
    // step 4
    Complex_scalar_field phi_hat2(G_hat2.get_points().get_shape(),
				  G_hat2.get_physical_size(),
				  G_hat2.get_physical_offset());
    Int3 shape(G_hat2.get_points().get_shape());
    Double3 h(rho.get_cell_size());
    timer("misc");
    for(int i=0; i<G_hat2.get_points().get_length(); ++i) {
      phi_hat2.get_points().get_base_address()[i] = 	
	rho_hat2.get_points().get_base_address()[i]*
	G_hat2.get_points().get_base_address()[i]*
	h[0]*h[1]*h[2];
    }
    timer("calc phi_hat");
    return phi_hat2;
  }

  Real_scalar_field
  get_phi2(Real_scalar_field &rho, Complex_scalar_field &phi_hat2,
	   fftw_plan inv_plan)
  {
    // step 5
    Int3 num_points2(rho.get_points().get_shape());
    num_points2.scale(2);
    Real_scalar_field phi2(num_points2.vector(),
			   phi_hat2.get_physical_size(),
			   phi_hat2.get_physical_offset());
    timer("misc");
    fftw_execute_dft_c2r(inv_plan,
			reinterpret_cast<fftw_complex *>
			(phi_hat2.get_points().get_base_address()),
			phi2.get_points().get_base_address());
    timer("invfft phi");
    double norm = 1.0/phi2.get_points().get_length();
    phi2.get_points().scale(norm);
    timer("calc norm");
    return phi2;
  }

  Real_scalar_field
  get_phi(Real_scalar_field &rho, Real_scalar_field &phi2)
  {
    // step 6
    Real_scalar_field phi(rho.get_points().get_shape(),rho.get_physical_size(),
			  rho.get_physical_offset());
    Int3 shape(phi.get_points().get_shape());
    Int3 point;
    timer("misc");
    for(int i=0; i<shape[0]; ++i) {
      point[0] = i;
      for(int j=0; j<shape[1]; ++j) {
	point[1] = j;
	for(int k=0; k<shape[2]; ++k) {
	  point[2] = k;
	  phi.get_points().set(point,phi2.get_points().get(point));
	}
      }
    }
    timer("calc phi");
    return phi;	
  }

  void
  get_plans(Real_scalar_field &rho, fftw_plan &plan, fftw_plan &inv_plan)
  {
    Int3 shape(rho.get_points().get_shape());
    shape.scale(2);
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

  Real_scalar_field
  solver_fftw3_open(Real_scalar_field &rho, bool z_periodic)
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

    double t0 = time();
    reset_timer();
    fftw_plan plan, inv_plan;
    get_plans(rho,plan,inv_plan);
    timer("get plans");
    Complex_scalar_field rho_hat2 = get_rho_hat2(rho,plan);
    Complex_scalar_field G_hat2 = get_G_hat2(rho,z_periodic,plan);
    Complex_scalar_field phi_hat2 = get_phi_hat2(rho, rho_hat2,
						 G_hat2);
    timer("misc");
    Real_scalar_field phi2 = get_phi2(rho,phi_hat2,inv_plan);
    timer("misc");
    Real_scalar_field phi = get_phi(rho,phi2);
    timer("misc");
    std::cout << "time total: " << time() - t0 << std::endl;
    return phi;
  }

  double
  fft_tester(int nx, int ny, int nz)
  {
    int array_length = nx*ny*nz;
    int dims[3] = {nx,ny,nz};

    Double3 one(1.0,1.0,1.0);
    Double3 zero(0.0,0.0,0.0);
    Real_scalar_field v(dims,one,zero);
    Complex_scalar_field vhat(Int3(nx,ny,nz),one,zero);
    Real_scalar_field vhathat(dims,one,zero);

    Int3 index;
    for (index[0] = 0; index[0]<nx; ++index[0]) {
      for (index[1] = 0; index[1]<ny; ++index[1]) {
	for (index[2] = 0; index[2]<nz; ++index[2]) {
	  v.get_points().set(index,index[0]*100.0+index[1]*10+index[0]);
	}
      }
    }
    v.get_points().print("v");
//     double a[nx*ny*nz];
//     double ahathat[nx*ny*nz];
//     fftw_complex ahat[nx*ny*(nz/2+1)];

    double *a, *ahathat;
    fftw_complex *ahat;
    a = (double *)fftw_malloc(nx*ny*nz*sizeof(double));
    ahathat = (double *)fftw_malloc(nx*ny*nz*sizeof(double));
    ahat = (fftw_complex *)fftw_malloc(nx*ny*(nz/2+1)*sizeof(fftw_complex));
    fftw_plan plan, inv_plan;
    for (int i=0; i<nx*ny*nz; ++i) {
      a[i] = i+17.0;
    }
    plan = fftw_plan_dft_r2c_3d(nx,ny,nz,
			     a,
			     ahat,
			     FFTW_MEASURE);
//     plan = fftw_plan_dft_r2c(3,
// 			     dims,
// 			     v.get_points().get_base_address(),
// 			     reinterpret_cast<fftw_complex *>
// 			     (vhat.get_points().get_base_address()),
// 			     FFTW_MEASURE);
//    v.get_points().print("v");    
    std::cout << "about to execute " << plan << "\n";
    fftw_execute(plan);
    
//     double max_err(0.0);
//     double ratio;
//     for(int i=0; i<dims[0]; ++i) {
//       for(int j=0; j<dims[1]; ++j) {
// 	for(int k=0; k<dims[2]; ++k) {
// 	  if (a2[i][j][k] != 0.0) {
// 	    ratio = a[i][j][k].real()/a2[i][j][k].real();
// 	  } else {
// 	    ratio = a[i][j][k].real() + 1.0;
// 	  }
// 	  double abs_diff = fabs(ratio-1.0);
// 	  if (abs_diff > max_err) {
// 	    max_err = abs_diff;
// 	  }
// 	}
//       }
//     }
//     return max_err;
    return -1.0;
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
      try {
	kick = tau * factor * E.get_val(Double3(mbs.local_particles(0,n),
						mbs.local_particles(2,n),
						mbs.local_particles(4,n)));
      } catch(std::out_of_range e) {
	//       std::cout << "particle " << n << " out of range ("
	// 		<< mbs.local_particles(0,n) << ", "
	// 		<< mbs.local_particles(2,n) << ", "
	// 		<< mbs.local_particles(4,n) << ") "
	// 		<<"\n";
	kick = 0.0;
      }
      mbs.local_particles(index,n) -= kick;
    }
  }

  void
  full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs)
  {
    for(int axis=0; axis<3; ++axis) {
      Real_scalar_field E = solver_fftw3::calculate_E_n(phi,axis);
      solver_fftw3::apply_E_n_kick(E,axis,tau,mbs);
    }
  }
}
