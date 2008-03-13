#include "cylindrical.h"
#include "field_domain.h"
#include "array_nd/array_3d.h"
#include "math_constants.h"

void 
get_cylindrical_coords(Macro_bunch_store &mbs, Array_2d<double> &coords)
{
    for (int n = 0; n < mbs.local_num; ++n) {
        double x = mbs.local_particles(0,n);
        double y = mbs.local_particles(2,n);        
        double r = sqrt(x*x + y*y);
        double theta;
        if (x>=0.0) {
            if(y>=0.0) {
                theta = asin(y/r);
            } else {
                theta = 2*pi + asin(y/r);
            }
        } else {
            theta = pi - asin(y/r);
        }
        coords(0,n) = r;
        coords(1,n) = theta;
        coords(2,n) = mbs.local_particles(4,n); // z
    }
}

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
void
deposit_charge_cic_cylindrical(const Field_domain &fdomain, 
    Array_3d<double > &rho , const Macro_bunch_store& mbs,
    const Array_2d<double> &coords)
{
    std::vector<int> indices(3);
    std::vector<double> offsets(3);
    std::vector<double> cell_size(fdomain.get_cell_size());
    std::vector<bool> periodic(fdomain.get_periodic());
    std::vector<int> grid_shape(fdomain.get_grid_shape());
    rho.set_all(0.0);
    for (int n = 0; n < mbs.local_num; ++n) {
        double r = coords(0,n);
        double theta = coords(1,n);
        double z = coords(2,n);
        fdomain.get_leftmost_indices_offsets(r,theta,z,indices,offsets);
        for (int i = 0; i < 2; ++i) {
            int this_i = indices[0] + i;
            if (periodic[0]) {
                this_i = this_i % grid_shape[0];
                if (this_i < 0) {
                    this_i += grid_shape[0];
                }
            }
            for (int j = 0; j < 2; ++j) {
                int this_j = indices[1] + j;
                if (periodic[1]) {
                    this_j = this_j % grid_shape[1];
                    if (this_j < 0) {
                        this_j += grid_shape[1];
                    }
                }
                for (int k = 0; k < 2; ++k) {
                    int this_k = indices[2] + k;
                    if (periodic[2]) {
                        this_k = this_k % grid_shape[2];
                        if (this_k < 0) {
                            this_k += grid_shape[2];
                        }
                    }
                    double weight = (1 - i - (1 - 2 * i) * offsets[0]) *
                                    (1 - j - (1 - 2 * j) * offsets[1]) *
                                    (1 - k - (1 - 2 * k) * offsets[2]);
        
                    if (rho.bounds_check(this_i,this_j,this_k)) {
                        rho(this_i,this_j,this_k) += weight;
                    }
                }
            }
        }
    }
}

