#include <mpi.h>
#include "deposit.h"
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>

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
    Int3 shape(sf.get_points().get_shape());
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
                 if(((indices[0] + i) != 0) && ((indices[0] + i) != shape[0]-1) // on the transverse grid edges 
                    && ((indices[1] + j) != 0) && ((indices[1] + j) != shape[1]-1)){   // no charge is deposited 
		     sf.get_points().add_to_point(Int3(indices[0] + i,
                                                      indices[1] + j,
                                                      indices[2] + k),
                                                 weight);
                    total_charge_per_cell_vol += weight;
                  }
                }
            }
        }
    }
 

      
       
     for (int i = 1; i < shape[0]-1; ++i) {
       for (int j = 1; j < shape[1]-1; ++j) {
            Int3 left(i, j, 0), right(i, j, shape[2] - 1);
            double sum;
	    sum = sf.get_points().get(left) +
                      sf.get_points().get(right); 
            if (!z_periodic) {                             
		total_charge_per_cell_vol -= sum;
		sum=0.;}	      
              sf.get_points().set(left, sum);
              sf.get_points().set(right, sum);
       }
    }
 

   //  std::cout<<"total charge is "<<total_charge_per_cell_vol * h[0]*h[1]*h[2]<<" cel size="<<h[0]*h[1]*h[2]<<std::endl<<std::endl;
//          }
    //jfa next line is wrong on multiple processors!!!
    // AM:  the total charge deposit on the grid is not equal to the number of particles 
    //when are particles outside the grid...	
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



