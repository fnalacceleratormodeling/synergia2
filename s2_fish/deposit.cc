#include "deposit.h"
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <mpi.h>

#include <cstdlib>

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
double
deposit_charge_cic(Real_scalar_field& sf, Macro_bunch_store& mbs,
                   bool z_periodic)
{
    Int3 indices;
    Double3 offsets;
    double total_charge_per_cell_vol = 0.0;
    Double3 h(sf.get_cell_size());
    double weight0 = 1.0 / (h[0] * h[1] * h[2]);
    //~ weight0 = 1.0; //jfa DEBUG ONLY!!!!!!!!!!!!!!!!!!!!!
    //~ std::cout << "jfa: btw, zperiodic = " << z_periodic << std::endl;
    sf.get_points().zero_all();
    for (int n = 0; n < mbs.local_num; ++n) {
        Double3 location(mbs.local_particles(0, n),
                         mbs.local_particles(2, n),
                         mbs.local_particles(4, n));
        indices = sf.get_leftmost_indices(location);
        offsets = sf.get_leftmost_offsets(location);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    double weight = weight0 * (1 - i - (1 - 2 * i) * offsets[0]) *
                                    (1 - j - (1 - 2 * j) * offsets[1]) *
                                    (1 - k - (1 - 2 * k) * offsets[2]);
                    sf.get_points().add_to_point(Int3(indices[0] + i,
                                                      indices[1] + j,
                                                      indices[2] + k),
                                                 weight);
                    total_charge_per_cell_vol += weight;
                }
            }
        }
    }
    if (z_periodic) {
        Int3 shape(sf.get_points().get_shape());
        double sum;
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                Int3 left(i, j, 0), right(i, j, shape[2] - 1);
                sum = sf.get_points().get(left) +
                      sf.get_points().get(right);
                sf.get_points().set(left, sum);
                sf.get_points().set(right, sum);				
            }
        }
    }

// 	Int3 shape(sf.get_points().get_shape());
// 	for (int j = 0; j < shape[1]; ++j) {
// 	   for (int k = 0; k < shape[2]-1; ++k) {
// 	Int3 left(0, j, k), right(shape[0]-1,j,k);
//         Int3 left1(1, j, k), right1(shape[0]-2,j,k);
// 	if((sf.get_points().get(left)>1.e-10)||(sf.get_points().get(right)>1.e-10))
//         std::cout<<"charge left="<<sf.get_points().get(left)<<"charge right="<<sf.get_points().get(right)<<"   jk= "<<j<<k<<std::endl;
// //           sf.get_points().set(right, 0.);
// //  	  sf.get_points().set(left, 0.); 
// // 	  sf.get_points().set(right1, 0.);
// //  	  sf.get_points().set(left1, 0.);
// 
//            }
//          }
    //jfa next line is wrong on multiple processors!!!
    return total_charge_per_cell_vol * h[0]*h[1]*h[2];
}

// Deposit charge on the nearest grid point to the particle
double
deposit_charge_ngp(Real_scalar_field& sf, Macro_bunch_store& mbs)
{
    Int3 indices;
    Double3 offsets;
    double total_charge_per_cell_vol = 0.0;
    Double3 h(sf.get_cell_size());
    double weight = 1.0 / (h[0] * h[1] * h[2]);
    sf.get_points().zero_all();
    for (int n = 0; n < mbs.local_num; ++n) {
        Double3 location(mbs.local_particles(0, n),
                         mbs.local_particles(2, n),
                         mbs.local_particles(4, n));
        indices = sf.get_leftmost_indices(location);
        offsets = sf.get_leftmost_offsets(location);
        for (int i = 0; i < 3; ++i) {
            if (offsets[i] > 0.5) {
                indices[i] += 1;
            }
        }
        sf.get_points().add_to_point(indices, weight);
        total_charge_per_cell_vol += weight;
    }
    return total_charge_per_cell_vol * h[0]*h[1]*h[2];
}

void
rho_to_rwvars(Real_scalar_field &rho, Array_1d<double> &zdensity,
            Array_1d<double> &xmom, Array_1d<double> &ymom)
{
	// returns zdensity in units of (# macroparticles),
	// xmom and ymom in units of meters 
    Int3 shape = rho.get_points().get_shape();
    std::vector<double> cell_size = rho.get_cell_size();
    double cell_volume = cell_size[0]*cell_size[1]*cell_size[2];
    std::vector<double> left = rho.get_left();
    for (int k=0; k<shape[2]; ++k) {
        zdensity(k) = 0.0;
        xmom(k) = 0.0;
        ymom(k) = 0.0;
        for (int i=0; i<shape[0]; ++i) {
            double x = left[0] + i*cell_size[0];
            std::cout << "jfa x = " << x << std::endl;
            for (int j=0; j<shape[1]; ++j) {
                double y = left[1] + j*cell_size[1];
                zdensity(k) += cell_volume*rho.get_points().get(Int3(i,j,k));
                xmom(k) += x*cell_volume*rho.get_points().get(Int3(i,j,k));
                ymom(k) += y*cell_volume*rho.get_points().get(Int3(i,j,k));
            }
        }
    }
    for (int k=0; k<shape[2]; ++k) {
        if (zdensity(k) > 0.0) {
            xmom(k) *= 1.0/zdensity(k);
            ymom(k) *= 1.0/zdensity(k);
        } else {
            xmom(k) = 0.0;
            ymom(k) = 0.0;
        }
    }
}

void
calculate_rwvars(Macro_bunch_store& mbs,
                   Array_1d<double> &zdensity,
                   Array_1d<double> &xmom, Array_1d<double> &ymom,
                   double z_left, double z_length)
{
//	const double really_big = 1.0e30;
//	double zmin = really_big;
//	double zmax = -really_big;
//    for (int n = 0; n < mbs.local_num; ++n) {
//    	double z = mbs.local_particles(4,n);
//    	if (z < zmin) {
//    		zmin = z;
//    	} else if (z > zmax) {
//    		zmax = z;
//    	}
//	}
//    double z_left = zmin;
//    double z_length = (zmax - zmin);
//    std::cout << "jfa: z_left = " << z_left << ", z_length = " << z_length << std::endl;
	int z_num = zdensity.get_length();
    double h = z_length/z_num;
    Array_1d<double> local_zdensity(z_num);
    Array_1d<double> local_xmom(z_num);
    Array_1d<double> local_ymom(z_num);
    local_zdensity.set_all(0.0);
    local_xmom.set_all(0.0);
    local_ymom.set_all(0.0);
    for (int n = 0; n < mbs.local_num; ++n) {
        int bin = static_cast<int>((mbs.local_particles(4,n)-z_left)/h);
        if ((bin < z_num) && (bin >= 0)) {
            local_zdensity(bin) += 1;
            local_xmom(bin) += mbs.local_particles(0,n);
            local_ymom(bin) += mbs.local_particles(2,n);
        }
    }
    
    // jfa: the other deposit functions do not do the communication. This one does.
    // Obviously, something should change.
    MPI_Allreduce(reinterpret_cast<void*>(local_zdensity.get_data_ptr()),
                    reinterpret_cast<void*>(zdensity.get_data_ptr()),
                    z_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(reinterpret_cast<void*>(local_xmom.get_data_ptr()),
                    reinterpret_cast<void*>(xmom.get_data_ptr()),
                    z_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(reinterpret_cast<void*>(local_ymom.get_data_ptr()),
                    reinterpret_cast<void*>(ymom.get_data_ptr()),
                    z_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (int k = 0; k < z_num; ++k) {
    	if (zdensity(k) != 0.0) {
	        xmom(k) /= zdensity(k)*mbs.units(0);
	        ymom(k) /= zdensity(k)*mbs.units(2);
    	} else {
    		xmom(k) = 0.0;
    		ymom(k) = 0.0;
    	}
    }
}   
