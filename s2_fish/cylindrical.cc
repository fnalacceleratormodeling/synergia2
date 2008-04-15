#include "cylindrical.h"
#include "math_constants.h"
#include <fftw3.h>

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

inline int
get_this(int index, int which, const std::vector<int> &indices, 
    const std::vector<int> &grid_shape, const std::vector<bool> &periodic)
{
    int this_ = indices[which] + index;
    if (periodic[which]) {
        this_ = this_ % grid_shape[which];
        if (this_ < 0) {
            this_ += grid_shape[which];
        }
    }
    return this_;
}

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
void
deposit_charge_cic_cylindrical(const Field_domain &fdomain, 
    Array_3d<double > &rho , Macro_bunch_store& mbs,
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
        for (int i = 0; i < 1; ++i) {
            int this_i = get_this(i,0,indices,grid_shape,periodic);
            for (int j = 0; j < 1; ++j) {
                int this_j = get_this(j,1,indices,grid_shape,periodic);
                for (int k = 0; k < 1; ++k) {
                    int this_k = get_this(k,2,indices,grid_shape,periodic);
                    double weight = 
                                    (1 - i - (1 - 2 * i) * offsets[0]*offsets[0]) *
                                    (1 - j - (1 - 2 * j) * offsets[1]) *
                                    (1 - k - (1 - 2 * k) * offsets[2]);
                    // jfa: The good news: the next line is correct.
                    // jfa: The bad news: I don't understand why.
                    // jfa: The 0.75 is a little bit of a mystery.
                    weight *= static_cast<double>(grid_shape[0])/
                        (2*indices[0]+0.75); // Jacobian factor
                    if (rho.bounds_check(this_i,this_j,this_k)) {
                        rho(this_i,this_j,this_k) += weight;
                    } else {
                        std::cout << "jfa debug oob: " << this_i << ", " << this_j << ", " << this_k << std::endl;
                    }
                }
            }
        }
    }
    if (periodic[0]) {
        for(int j=0; j<grid_shape[1]; ++j) {
            for(int k=0; k<grid_shape[2]; ++k) {
                double sum = rho(0,j,k) + rho(grid_shape[0]-1,j,k);
                rho(0,j,k) = sum;
                rho(grid_shape[0]-1,j,k) = sum;
            }
        }
    }
    if (periodic[1]) {
        for(int i=0; i<grid_shape[0]; ++i) {
            for(int k=0; k<grid_shape[2]; ++k) {
                double sum = rho(i,0,k) + rho(i,grid_shape[1]-1,k);
                rho(i,0,k) = sum;
                rho(i,grid_shape[1]-1,k) = sum;
            }
        }
    }
    if (periodic[2]) {
        for(int i=0; i<grid_shape[0]; ++i) {
            for(int j=0; j<grid_shape[1]; ++j) {
                double sum = rho(i,j,0) + rho(i,j,grid_shape[2]-1);
                rho(i,j,0) = sum;
                rho(i,j,grid_shape[2]-1) = sum;
            }
        }
    }
}

// Adapted from a GSL routine.
// plain gauss elimination, only not bothering with the zeroes
//
//       diag[0]  abovediag[0]             0   .....
//  belowdiag[0]       diag[1]  abovediag[1]   .....
//             0  belowdiag[1]       diag[2]
//             0             0  belowdiag[2]   .....
//
void
solve_tridiag_nonsym(const Array_1d<double> &diag,
    const Array_1d<double> &abovediag,
    const Array_1d<double> &belowdiag,
    const Array_1d<double> &rhs,
    Array_1d<double> &x)
{
    int N = diag.get_shape()[0];
    Array_1d<double> alpha(N);
    Array_1d<double> z(N);
    size_t i, j;

    // Bidiagonalization (eliminating belowdiag)
    // & rhs update
    // diag' = alpha
    // rhs' = z
    alpha(0) = diag(0);
    z(0) = rhs(0);

    for (i = 1; i < N; ++i) {
        const double t = belowdiag(i - 1)/alpha(i-1);
        alpha(i) = diag(i) - t*abovediag(i - 1);
        z(i) = rhs(i) - t*z(i-1);
        if (alpha(i) == 0) {
            throw
                std::runtime_error("solve_tridiag_nonsym: zero pivot");
        }
    }

    // backsubstitution
    x(N - 1) = z(N - 1)/alpha(N - 1);
    if (N >= 2) {
      for (i = N - 2, j = 0; j <= N - 2; ++j, --i) {
          x(i) = (z(i) - abovediag(i) * x(i + 1))/alpha(i);
        }
    }
}

void 
solve_cylindrical_finite_periodic(const Field_domain &fdomain,
    Array_3d<double > &rho, Array_3d<double> &phi)
{
    std::cout << "\njfa wtf:" << fdomain.get_grid_shape()[0] << std::endl;
    std::cout << rho.get_strides()[0] << " " << rho.get_strides()[1]
        << " " << rho.get_strides()[2] << std::endl;
    std::vector<int> shape = fdomain.get_grid_shape();
    Array_3d<std::complex<double> > rho_lm(
        vector3(shape[0],shape[1],shape[2]/2+1));
        
    fftw_plan plan = fftw_plan_many_dft_r2c(2,
        &shape[1],
        shape[0],
        rho.get_data_ptr(),
        NULL,
        1,
        shape[0],
        reinterpret_cast<double (*)[2]>(rho_lm.get_data_ptr()),
        NULL,
        1,
        shape[0],
        FFTW_ESTIMATE);
    fftw_execute(plan);
    
 }
