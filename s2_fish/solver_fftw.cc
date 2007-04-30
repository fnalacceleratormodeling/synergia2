#include "solver_fftw.h"
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

#undef DL_IMPORT
#include "mytimer.h"
#include "fftw2_helper.h"

#include "communicate.h"

Complex_scalar_field
get_rho_hat2(Real_scalar_field &rho, Fftw_helper &fftwh) {
	// steps 1 and 2
	Int3 num_points2 = rho.get_points().get_shape();
	num_points2.scale(2);
	Double3 physical_size2 = rho.get_physical_size();
	physical_size2.scale(2.0);
	Real_scalar_field rho2(fftwh.padded_shape_real().vector(),
	                       physical_size2.vector(),
	                       rho.get_physical_offset(),
	                       fftwh.guard_lower(), fftwh.guard_upper());
	//  rho2.get_points().set_storage_size(fftwh.local_size());
	timer("misc");
	rho2.get_points().zero_all();
	timer("calc zero all rho2");
	Int3 num_points = rho.get_points().get_shape();
	Int3 index;
	timer("misc");
	int index0_max = std::min(fftwh.upper(), num_points[0]);
	for (index[0] = fftwh.lower(); index[0] < index0_max; ++index[0]) {
		for (index[1] = 0; index[1] < num_points[1]; ++index[1]) {
			for (index[2] = 0; index[2] < num_points[2]; ++index[2]) {
				rho2.get_points().set(index, rho.get_points().get(index));
			}
		}
	}
	timer("calc rho2");

	Complex_scalar_field rho_hat2(fftwh.padded_shape_complex().vector(),
	                              physical_size2.vector(),
	                              rho.get_physical_offset(),
	                              fftwh.guard_lower(), fftwh.guard_upper());
	//  rho_hat2.get_points().set_storage_size(fftwh.local_size());
	timer("rho plan");
	fftwh.transform(rho2, rho_hat2);
	timer("rho fft");
	return rho_hat2;
}

Real_scalar_field
get_G2(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh) {
	const double pi = 4.0 * atan(1.0);
	Int3 num_points = rho.get_points().get_shape();
	Int3 num_points2 = rho.get_points().get_shape();
	num_points2.scale(2);
	Double3 physical_size = rho.get_physical_size();
	Double3 physical_size2 = rho.get_physical_size();
	physical_size2.scale(2.0);
	Real_scalar_field G2(fftwh.padded_shape_real().vector(),
	                     physical_size2.vector(),
	                     rho.get_physical_offset(),
	                     fftwh.guard_lower(), fftwh.guard_upper());
	//  G2.get_points().set_storage_size(fftwh.local_size());
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
	double G000 = (1.0 / 4.0 * pi) * (3.0 / (2.0 * (sqrt(3)) *
	                                  sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2])));
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
	double x, y, z, G;
	const int num_images = 4;
	int miy, miz; // mirror index x, etc.
	for (index[0] = fftwh.lower(); index[0] < fftwh.upper(); ++index[0]) {
		if (index[0] > num_points2[0] / 2) {
			x = (num_points2[0] - index[0]) * h[0];
		} else {
			x = index[0] * h[0];
		}
		for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
			y = index[1] * h[1];
			miy = num_points2[1] - index[1];
			if (miy == num_points[1]) {
				miy = num_points2[1];
			}
			for (index[2] = 0; index[2] <= num_points[2]; ++index[2]) {
				z = index[2] * h[2];
				miz = num_points2[2] - index[2];
				if (miz == num_points[2]) {
					miz = num_points2[2];
				}
				if (!((x == 0.0) && (y == 0.0) && (z == 0.0))) {
					G = 1.0 / (4.0 * pi * sqrt(x * x + y * y + z * z));
				} else {
					G = G000;
				}
				if (z_periodic) {
					for (int image = -num_images; image <= num_images; ++image) {
						if (image != 0) {
							double z_image = z + image * physical_size[2];
							if (!((x == 0.0) && (y == 0.0) && (fabs(z_image) < 1.0e-14))) {
								G += 1.0 / (4.0 * pi * sqrt(x * x + y * y + z_image * z_image));
							} else {
								G += G000;
							}
						}
					}
				}
				G2.get_points().set(index, G);
				// three mirror images
				if (miy < num_points2[1]) {
					G2.get_points().set(Int3(index[0], miy, index[2]), G);
					if (miz < num_points2[2]) {
						G2.get_points().set(Int3(index[0], miy, miz), G);
					}
				}
				if (miz < num_points2[2]) {
					G2.get_points().set(Int3(index[0], index[1], miz), G);
				}
			}
		}
	}
	timer("calc G");
	return G2;
}

Real_scalar_field
get_G2_old(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh) {
	const double pi = 4.0 * atan(1.0);
	Int3 num_points = rho.get_points().get_shape();
	Int3 num_points2 = rho.get_points().get_shape();
	num_points2.scale(2);
	Double3 physical_size = rho.get_physical_size();
	Double3 physical_size2 = rho.get_physical_size();
	physical_size2.scale(2.0);
	Real_scalar_field G2(fftwh.padded_shape_real().vector(),
	                     physical_size2.vector(),
	                     rho.get_physical_offset(),
	                     fftwh.guard_lower(), fftwh.guard_upper());
	//  G2.get_points().set_storage_size(fftwh.local_size());
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
	double G000 = (1.0 / 4.0 * pi) * (3.0 / (2.0 * (sqrt(3)) *
	                                  sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2])));
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
	double x, y, z, G;
	const int num_images = 4;
	int mix, miy, miz; // mirror index x, etc.
	for (index[0] = 0; index[0] <= num_points[0]; ++index[0]) {
		x = index[0] * h[0];
		mix = num_points2[0] - index[0];
		if (mix == num_points[0]) {
			mix = num_points2[0];
		}
		for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
			y = index[1] * h[1];
			miy = num_points2[1] - index[1];
			if (miy == num_points[1]) {
				miy = num_points2[1];
			}
			for (index[2] = 0; index[2] <= num_points[2]; ++index[2]) {
				z = index[2] * h[2];
				miz = num_points2[2] - index[2];
				if (miz == num_points[2]) {
					miz = num_points2[2];
				}
				if (!((x == 0.0) && (y == 0.0) && (z == 0.0))) {
					G = 1.0 / (4.0 * pi * sqrt(x * x + y * y + z * z));
				} else {
					G = G000;
				}
				if (z_periodic) {
					for (int image = -num_images; image <= num_images; ++image) {
						if (image != 0) {
							double z_image = z + image * physical_size[2];
							if (!((x == 0.0) && (y == 0.0) && (fabs(z_image) < 1.0e-14))) {
								G += 1.0 / (4.0 * pi * sqrt(x * x + y * y + z_image * z_image));
							} else {
								G += G000;
							}
						}
					}
				}
				G2.get_points().set(index, G);
				// seven mirror images
				if (mix < num_points2[0]) {
					G2.get_points().set(Int3(mix, index[1], index[2]), G);
					if (miy < num_points2[1]) {
						G2.get_points().set(Int3(mix, miy, index[2]), G);
						if (miz < num_points2[2]) {
							G2.get_points().set(Int3(mix, miy, miz), G);
						}
					}
					if (miz < num_points2[2]) {
						G2.get_points().set(Int3(mix, index[1], miz), G);
					}
				}
				if (miy < num_points2[1]) {
					G2.get_points().set(Int3(index[0], miy, index[2]), G);
					if (miz < num_points2[2]) {
						G2.get_points().set(Int3(index[0], miy, miz), G);
					}
				}
				if (miz < num_points2[2]) {
					G2.get_points().set(Int3(index[0], index[1], miz), G);
				}
			}
		}
	}
	timer("calc G");
	return G2;
}

Complex_scalar_field
get_G_hat2(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh) {
	//step 3
	Real_scalar_field G2 = get_G2(rho, z_periodic, fftwh);
	Complex_scalar_field G_hat2(fftwh.padded_shape_complex().vector(),
	                            G2.get_physical_size(),
	                            rho.get_physical_offset(),
	                            fftwh.guard_lower(), fftwh.guard_upper());
	//  G_hat2.get_points().set_storage_size(fftwh.local_size());
	fftwh.transform(G2, G_hat2);
	timer("G fft");
	return G_hat2;
}

Complex_scalar_field
get_phi_hat2(Real_scalar_field &rho, Complex_scalar_field &rho_hat2,
             Complex_scalar_field &G_hat2, Fftw_helper &fftwh) {
	// step 4
	Complex_scalar_field phi_hat2(fftwh.padded_shape_complex().vector(),
	                              G_hat2.get_physical_size(),
	                              G_hat2.get_physical_offset(),
	                              fftwh.guard_lower(), fftwh.guard_upper());
	//  phi_hat2.get_points().set_storage_size(fftwh.local_size());
	Int3 shape(G_hat2.get_points().get_shape());
	Double3 h(rho.get_cell_size());
	timer("misc");
	for (int i = 0; i < G_hat2.get_points().get_length(); ++i) {
		phi_hat2.get_points().get_base_address()[i] =
		    rho_hat2.get_points().get_base_address()[i] *
		    G_hat2.get_points().get_base_address()[i] *
		    h[0] * h[1] * h[2];
	}
	timer("calc phi_hat");
	return phi_hat2;
}

Real_scalar_field
get_phi2(Real_scalar_field &rho, Complex_scalar_field &phi_hat2,
         Fftw_helper &fftwh) {
	// step 5
	Int3 num_points2(rho.get_points().get_shape());
	num_points2.scale(2);
	Real_scalar_field phi2(fftwh.padded_shape_real().vector(),
	                       phi_hat2.get_physical_size(),
	                       phi_hat2.get_physical_offset(),
	                       fftwh.guard_lower(), fftwh.guard_upper());
	//  phi2.get_points().set_storage_size(fftwh.local_size());
	timer("misc");
	fftwh.inv_transform(phi_hat2, phi2);
	timer("invfft phi");
	double norm = 1.0 / (num_points2[0] * num_points2[1] * num_points2[2]);
	phi2.get_points().scale(norm);
	timer("calc norm");
	return phi2;
}

Real_scalar_field
get_phi(Real_scalar_field &rho, Real_scalar_field &phi2, Fftw_helper &fftwh) {
	// step 6
	Real_scalar_field phi(rho.get_points().get_shape(), rho.get_physical_size(),
	                      rho.get_physical_offset(),
	                      fftwh.guard_lower(),
	                      std::min(fftwh.guard_upper(),
	                               rho.get_points().get_shape()[0]));
	Int3 shape(phi.get_points().get_shape());
	Int3 point;
	timer("misc");
	int i_max = std::min(fftwh.upper(), shape[0]);
	for (int i = fftwh.lower(); i < i_max; ++i) {
		point[0] = i;
		for (int j = 0; j < shape[1]; ++j) {
			point[1] = j;
			for (int k = 0; k < shape[2]; ++k) {
				point[2] = k;
				phi.get_points().set(point, phi2.get_points().get(point));
			}
		}
	}
	timer("calc phi");
	return phi;
}

void
fill_guards(Real_scalar_field &rho, Fftw_helper &fftwh)
{
    Int3 shape(rho.get_points().get_shape());
    size_t message_size = shape[1]*shape[2];
    void *recv_buffer, *send_buffer;
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;

    // send to the right
    if(fftwh.lower() == fftwh.guard_lower()) {
        recv_buffer = malloc(message_size*sizeof(double));
    } else {
        recv_buffer = reinterpret_cast<void*>(rho.get_points().get_offset_base_address(fftwh.guard_lower()));
    }
    send_buffer = reinterpret_cast<void*>(rho.get_points().get_offset_base_address(fftwh.upper()-1));
    if (rank < size -1) {
            MPI_Send(send_buffer,message_size,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD);
    }
    if (rank > 0) {
            MPI_Recv(recv_buffer,message_size,MPI_DOUBLE,rank-1,rank-1,MPI_COMM_WORLD,&status);
    }
    std::cout << "completed right on rank " << rank << std::endl;
    if(fftwh.lower() == fftwh.guard_lower()) {
        free(recv_buffer);
    }
    
    //send to the left
    if(fftwh.upper() == fftwh.guard_upper()) {
        recv_buffer = malloc(message_size*sizeof(double));
    } else {
        recv_buffer = reinterpret_cast<void*>(rho.get_points().get_offset_base_address(fftwh.guard_upper()-1));
    }
    send_buffer = reinterpret_cast<void*>(rho.get_points().get_offset_base_address(fftwh.lower()));
    if (rank > 0) {
            MPI_Send(send_buffer,message_size,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD);
    }
    if (rank < size - 1) {
            MPI_Recv(recv_buffer,message_size,MPI_DOUBLE,rank+1,rank+1,MPI_COMM_WORLD,&status);
    }
    std::cout << "completed left on rank " << rank << std::endl;
    if(fftwh.upper() == fftwh.guard_upper()) {
        free(recv_buffer);
    }
    }

Real_scalar_field
solver_fftw_open(Real_scalar_field &rho, bool z_periodic) {
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
    gather_rho(rho);
    timer("gather rho");
	Fftw2_helper_mpi fftwh(rho);
    timer("get plans");
	Complex_scalar_field rho_hat2 = get_rho_hat2(rho, fftwh);
	Complex_scalar_field G_hat2 = get_G_hat2(rho, z_periodic, fftwh);
	Complex_scalar_field phi_hat2 = get_phi_hat2(rho, rho_hat2, G_hat2, fftwh);
	timer("misc");
	Real_scalar_field phi2 = get_phi2(rho, phi_hat2, fftwh);
	timer("misc");
	Real_scalar_field phi = get_phi(rho, phi2, fftwh);
	timer("misc");
    fill_guards(rho,fftwh);
//   std::cout << "time total: " << time() - t0 << std::endl;
	return phi;
}

