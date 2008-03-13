#include "cylindrical.h"
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
//~ double
//~ deposit_charge_cic_cylindrical(Real_scalar_field& sf, Macro_bunch_store& mbs,
                   //~ bool z_periodic)
//~ {
    //~ Int3 indices;
    //~ Double3 offsets;
    //~ double total_charge_per_cell_vol = 0.0;
    //~ Double3 h(sf.get_cell_size());
    //~ double weight0 = 1.0 / (h[0] * h[1] * h[2]);
    //~ sf.get_points().zero_all();
    //~ for (int n = 0; n < mbs.local_num; ++n) {
        //~ Double3 location(mbs.local_particles(0, n),
                         //~ mbs.local_particles(2, n),
                         //~ mbs.local_particles(4, n));
        //~ indices = sf.get_leftmost_indices(location);
        //~ offsets = sf.get_leftmost_offsets(location);
        //~ for (int i = 0; i < 2; ++i) {
            //~ for (int j = 0; j < 2; ++j) {
                //~ for (int k = 0; k < 2; ++k) {
                    //~ double weight = weight0 * (1 - i - (1 - 2 * i) * offsets[0]) *
                                    //~ (1 - j - (1 - 2 * j) * offsets[1]) *
                                    //~ (1 - k - (1 - 2 * k) * offsets[2]);
                        //~ sf.get_points().add_to_point(Int3(indices[0] + i,
                                                          //~ indices[1] + j,
                                                          //~ indices[2] + k),
                                                     //~ weight);
                        //~ total_charge_per_cell_vol += weight;
                //~ }
            //~ }
        //~ }
    //~ }
    //~ if (z_periodic) {
        //~ Int3 shape(sf.get_points().get_shape());
        //~ double sum;
        //~ for (int i = 0; i < shape[0]; ++i) {
            //~ for (int j = 0; j < shape[1]; ++j) {
                //~ Int3 left(i, j, 0), right(i, j, shape[2] - 1);
                //~ sum = sf.get_points().get(left) +
                      //~ sf.get_points().get(right);
                //~ sf.get_points().set(left, sum);
                //~ sf.get_points().set(right, sum);
            //~ }
        //~ }
    //~ }
    //~ //jfa next line is wrong on multiple processors!!!
    //~ return total_charge_per_cell_vol * h[0]*h[1]*h[2];
//~ }

