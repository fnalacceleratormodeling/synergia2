#include "cylindrical.h"
#include "math_constants.h"
#include <fftw3.h>
#include <cmath>

void 
get_cylindrical_coords(Macro_bunch_store &mbs, Array_2d<double> &coords)
{
    for (int n = 0; n < mbs.local_num; ++n) {
        double x = mbs.local_particles(0,n);
        double y = mbs.local_particles(2,n);
        double r = sqrt(x*x + y*y);
        double theta;
        if (r == 0.0) {
            theta = 0.0;
        } else {
            if (x>=0.0) {
                if(y>=0.0) {
                    theta = asin(y/r);
                } else {
                    theta = 2*pi + asin(y/r);
                }
            } else {
                theta = pi - asin(y/r);
            }
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

int numout, numin;

void
add_to_cylindrical_cell(Array_3d<double > &rho,
                        int ir, int iphi, int iz,
                        double scaled_rmin, double scaled_rmax,
                        double scaled_overlap_phi, double scaled_overlap_z,
                        double cell_size_r, double cell_size_phi, double cell_size_z,
                        double cloud_volume)
{
    double overlap_volume = 0.5*
                (scaled_rmax*scaled_rmax - scaled_rmin*scaled_rmin)*
                scaled_overlap_phi*scaled_overlap_z*
                cell_size_r*cell_size_r*cell_size_phi*cell_size_z;
    double cell_volume = 0.5*((ir+1)*(ir+1) - ir*ir)*cell_size_r*cell_size_r*
                cell_size_phi*cell_size_z;
//         std::cout
//             << ", " << overlap_volume
//             << ", " << cell_volume
//             << ", " << cloud_volume
//             << ", " << overlap_volume/(cell_volume*cloud_volume) << std::endl;

    if (rho.bounds_check(ir,iphi,iz)) {
        rho(ir,iphi,iz) += overlap_volume/(cell_volume*cloud_volume);
        ++numin;
/*        std::cout << "add_to_cylindrical_cell: "
                << ir
                << ", " << iphi
                << ", " << iz 
                << " " <<  scaled_rmin
                << " " <<  scaled_rmax
                << " " <<  scaled_overlap_phi
                << " " <<  scaled_overlap_z
                << " " <<  overlap_volume/(cell_volume*cloud_volume)
                << std::endl;*/
    } else {
        std::cout << "add_to_cylindrical_cell outside bounds: " << ir
                << ", " << iphi
                << ", " << iz << std::endl;
        ++numout;
    }
}

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
void
deposit_charge_cic_cylindrical(const Cylindrical_field_domain &fdomain, 
    Array_3d<double > &rho, Macro_bunch_store& mbs,
    const Array_2d<double> &coords)
{
    std::vector<int> indices(3);
    std::vector<double> offsets(3);
// jfa: we assume z to be periodic!!!!!!!!    std::vector<bool> periodic(fdomain.get_periodic());
    std::vector<int> grid_shape(fdomain.get_grid_shape());
    std::vector<double> cell_size(fdomain.get_cell_size());
    
    rho.set_all(0.0);
    numout = 0;
    numin = 0;
    for (int n = 0; n < mbs.local_num; ++n) {
        double r = coords(0,n);
        double theta = coords(1,n);
        double z = coords(2,n);
        fdomain.get_leftmost_indices_offsets(r,theta,z,indices,offsets);
        
        int left_ir, left_iphi, left_iz;
        int right_ir,right_iphi, right_iz;
        double left_r, center_r, right_r;
        double left_overlap_phi, left_overlap_z;
        
        if (offsets[0] > 0.5) {
            left_ir = indices[0];
            right_ir = indices[0] + 1;
            left_r = indices[0] + offsets[0] - 0.5;
            center_r = indices[0] + 1;
            right_r = left_r + 1;
        } else {
            left_ir = indices[0] - 1;
            right_ir = indices[0];
            left_r = indices[0] + offsets[0] - 0.5;
            center_r = indices[0];
            right_r = left_r + 1;
        }
        // The cloud volume is the volume of the cell containing the particle
        double cloud_volume;
        if (left_ir < 0) {
            // use smallest cell volume
            cloud_volume = 0.5*cell_size[0]*cell_size[0]*cell_size[1]*cell_size[2];
        } else {
            // use volume of a cell centered on particle
            cloud_volume = 0.5*((indices[0]+offsets[0]+0.5)*(indices[0]+offsets[0]+0.5) - 
                        (indices[0]+offsets[0]-0.5)*(indices[0]+offsets[0]-0.5))
                        *cell_size[0]*cell_size[0]*cell_size[1]*cell_size[2];
        }

        
        if (offsets[1] > 0.5) {
            left_iphi = indices[1];
            if (indices[1] == grid_shape[1] - 1) {
                right_iphi = 0;
            } else {
                right_iphi = indices[1] + 1;
            }
            left_overlap_phi = 1 - offsets[1];
        } else {
            if (indices[1] == 0) {
                left_iphi = grid_shape[1] - 1;
            } else {
                left_iphi = indices[1] - 1;
            }
            right_iphi = indices[1];
            left_overlap_phi = offsets[1];
        }
        
        if (offsets[2] > 0.5) {
            left_iz = indices[2];
            if (indices[2] == grid_shape[2] - 1) {
                right_iz = 0;
            } else {
                right_iz = indices[2] + 1;
            }
            left_overlap_z = 1 - offsets[2];
        } else {
            if (indices[2] == 0) {
                left_iz = grid_shape[2] - 1;
            } else {
                left_iz = indices[2] - 1;
            }
            right_iz = indices[2];
            left_overlap_z = offsets[2];
        }
        if (left_iz > 100) {std::cout << "left_iz = " << left_iz << std::endl;}
        if (right_iz > 100) {std::cout << "right_iz = " << right_iz << " " << indices[2] << " " << offsets[2] << " " << z << std::endl;}

        add_to_cylindrical_cell(rho,right_ir,left_iphi,left_iz,
                                center_r,right_r,
                                left_overlap_phi,left_overlap_z,
                                cell_size[0],cell_size[1],cell_size[2],
                                cloud_volume);
        add_to_cylindrical_cell(rho,right_ir,right_iphi,left_iz,
                                center_r,right_r,
                                1.0-left_overlap_phi,left_overlap_z,
                                cell_size[0],cell_size[1],cell_size[2],
                                cloud_volume);
        add_to_cylindrical_cell(rho,right_ir,left_iphi,right_iz,
                                center_r,right_r,
                                left_overlap_phi,1.0-left_overlap_z,
                                cell_size[0],cell_size[1],cell_size[2],
                                cloud_volume);
        add_to_cylindrical_cell(rho,right_ir,right_iphi,right_iz,
                                center_r,right_r,
                                1.0-left_overlap_phi,1.0-left_overlap_z,
                                cell_size[0],cell_size[1],cell_size[2],
                                cloud_volume);
        if (left_ir < 0) {
            for (int i_phi=0; i_phi<grid_shape[1]; ++i_phi) {
                add_to_cylindrical_cell(rho,0,i_phi,left_iz,
                                        right_r,1,
                                        1.0/grid_shape[1],left_overlap_z,
                                        cell_size[0],cell_size[1],cell_size[2],
                                        cloud_volume);
                add_to_cylindrical_cell(rho,0,i_phi,right_iz,
                                        right_r,1,
                                        1.0/grid_shape[1],1.0-left_overlap_z,
                                        cell_size[0],cell_size[1],cell_size[2],
                                        cloud_volume);
            }
/*            add_to_cylindrical_cell(rho,0,left_iphi,left_iz,
                                    right_r,1,
                                    left_overlap_phi,left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
            add_to_cylindrical_cell(rho,0,right_iphi,left_iz,
                                    right_r,1,
                                    1.0-left_overlap_phi,left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
            add_to_cylindrical_cell(rho,0,left_iphi,right_iz,
                                    right_r,1,
                                    left_overlap_phi,1.0-left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
            add_to_cylindrical_cell(rho,0,right_iphi,right_iz,
                                    right_r,1,
                                    1.0-left_overlap_phi,1.0-left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);*/
        } else {
            add_to_cylindrical_cell(rho,left_ir,left_iphi,left_iz,
                                    left_r,center_r,
                                    left_overlap_phi,left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
            add_to_cylindrical_cell(rho,left_ir,right_iphi,left_iz,
                                    left_r,center_r,
                                    1.0-left_overlap_phi,left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
            add_to_cylindrical_cell(rho,left_ir,left_iphi,right_iz,
                                    left_r,center_r,
                                    left_overlap_phi,1.0-left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
            add_to_cylindrical_cell(rho,left_ir,right_iphi,right_iz,
                                    left_r,center_r,
                                    1.0-left_overlap_phi,1.0-left_overlap_z,
                                    cell_size[0],cell_size[1],cell_size[2],
                                    cloud_volume);
        }
    }
    std::cout << "numin = " << numin << ", numout = " << numout << std::endl;
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
solve_tridiag_nonsym(const Array_1d<std::complex<double> > &diag,
    const Array_1d<std::complex<double> > &abovediag,
    const Array_1d<std::complex<double> > &belowdiag,
    const Array_1d<std::complex<double> > &rhs,
    Array_1d<std::complex<double> > &x)
{
    int N = diag.get_shape()[0];
    Array_1d<std::complex<double> > alpha(N);
    Array_1d<std::complex<double> > z(N);
    size_t i, j;

    // Bidiagonalization (eliminating belowdiag)
    // & rhs update
    // diag' = alpha
    // rhs' = z
    alpha(0) = diag(0);
    z(0) = rhs(0);

    for (i = 1; i < N; ++i) {
        const std::complex<double> t = belowdiag(i - 1)/alpha(i-1);
        alpha(i) = diag(i) - t*abovediag(i - 1);
        z(i) = rhs(i) - t*z(i-1);
        if (alpha(i) == 0.0) {
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

extern void
    array_2d_to_octave_file(const Array_2d<double> &array, const std::string filename);

void
array_2d_to_octave_file_real(const Array_2d<std::complex<double> > &array, const std::string filename)
{
    std::ofstream stream(filename.c_str());
    for (int i = 0; i < array.get_shape()[0]; ++i) {
        for (int j = 0; j < array.get_shape()[1]; ++j) {
            stream << std::setprecision(16) << array(i, j).real();
            if (j == array.get_shape()[1] - 1) {
                stream << std::endl;
            } else {
                stream << " ";
            }
        }
    }
    stream.close();
}

void
array_2d_to_octave_file_imag(const Array_2d<std::complex<double> > &array, const std::string filename)
{
    std::ofstream stream(filename.c_str());
    for (int i = 0; i < array.get_shape()[0]; ++i) {
        for (int j = 0; j < array.get_shape()[1]; ++j) {
            stream << std::setprecision(16) << array(i, j).imag();
            if (j == array.get_shape()[1] - 1) {
                stream << std::endl;
            } else {
                stream << " ";
            }
        }
    }
    stream.close();
}

void 
solve_cylindrical_finite_periodic(const Field_domain &fdomain,
    Array_3d<double > &rho, Array_3d<double> &phi)
{
    std::vector<int> shape = fdomain.get_grid_shape();
    double z_length = fdomain.get_physical_size()[2];
    // the shape of the FFT'd array (shape_lm) is halved in the third 
    // dimension because of the peculiar (but efficient) way FFTW does
    // real-to-complex transforms.
    std::vector<int> shape_lm = vector3(shape[0],shape[1],shape[2]/2+1);
    Array_3d<std::complex<double> > rho_lm(shape_lm);

    fftw_plan plan = fftw_plan_many_dft_r2c(2,
        &shape[1], shape[0],
        rho.get_data_ptr(),
        NULL, 1, shape[1]*shape[2],
        reinterpret_cast<double (*)[2]>(rho_lm.get_data_ptr()),
        NULL, 1, shape_lm[1]*shape_lm[2],
        FFTW_ESTIMATE);
    fftw_execute(plan);

    Array_1d<std::complex<double> > diag(shape_lm[0]);
    Array_1d<std::complex<double> > above_diag(shape_lm[0]-1);
    Array_1d<std::complex<double> > below_diag(shape_lm[0]-1);
    Array_3d<std::complex<double> > phi_lm(shape_lm);
    double deltar = fdomain.get_physical_size()[0]/
        (fdomain.get_grid_shape()[0]+0.5);

    for(int l=0; l<shape_lm[1]; ++l) {
        for(int m=0; m<shape_lm[2]; ++m) {
            for(int i=0; i<shape_lm[0]; ++i) {
                double r = (i+0.5)*deltar;
                diag(i) = -2.0*(1.0/(deltar*deltar));
                int wavenumber_l = (l+shape_lm[1]/2) % shape_lm[1] - 
                    shape_lm[1]/2;
                int wavenumber_m = (m+shape_lm[2]/2) % shape_lm[2] - 
                    shape_lm[2]/2;
                diag(i) += - wavenumber_l*wavenumber_l/(r*r) - 
                    pow(2*pi*wavenumber_m/z_length,2);
                if (i<(shape_lm[0]-1)) {
                    above_diag.at(i) = 1.0/(deltar*deltar) + 1.0/(2*deltar*r);
                }
                if (i>0) {
                    below_diag.at(i-1) = 1.0/(deltar*deltar) - 1.0/(2*deltar*r);
                }
            }
        Array_1d<std::complex<double> > rhs =
            rho_lm.slice(vector3(Range(),Range(l),Range(m)));
        Array_1d<std::complex<double> > x = 
            phi_lm.slice(vector3(Range(),Range(l),Range(m)));
         solve_tridiag_nonsym(diag,above_diag,below_diag,
            rhs,x);
        }
    }

    plan = fftw_plan_many_dft_c2r(2,
        &shape[1], shape_lm[0],
        reinterpret_cast<double (*)[2]>(phi_lm.get_data_ptr()),
        NULL, 1, shape_lm[1]*shape_lm[2],
        phi.get_data_ptr(),
        NULL, 1, shape[1]*shape[2],
        FFTW_ESTIMATE);
    fftw_execute(plan);
    // FFTW transforms are not normalized. We need to apply the normalization
    // manually.
    phi.scale(1.0/(shape[1]*shape[2]));
}

void
calculate_E_cylindrical(const Field_domain &fdomain,
                          Array_3d<double> &phi,
                          Array_3d<double> &Ex,
                          Array_3d<double> &Ey,
                          Array_3d<double> &Ez)
{
    std::vector<int> shape = phi.get_shape();
    std::vector<double> size = fdomain.get_physical_size();
    double theta_step = 2.0*pi/(shape[1]-1);
    double z_step = size[2]/(shape[2]-1);
    double ordinary_r_step = size[0]*2.0/(2*shape[0]+1);
    for(int i_r = 0; i_r<shape[0]; ++i_r) {
        int r_left = i_r-1;
        int r_right = i_r+1;
        double r_step = ordinary_r_step;
        double r = (i_r-0.5)*ordinary_r_step;
        if(i_r == 0) {
            r_left = 0;
            r_step *= 0.5;
        }
        if(i_r == shape[0]-1) {
            r_right = shape[0]-1;
            r_step = 0.5;
        }
        for(int i_theta = 0; i_theta<shape[1]; ++i_theta) {
            int theta_left = i_theta-1;
            int theta_right = i_theta+1;
            double theta = (i_theta-1)*theta_step;
            if (theta_left == -1) {
                theta_left = shape[1]-1;
            }
            if (theta_right == shape[1]) {
                theta_right = 0;
            }
            for(int i_z = 0; i_z<shape[2]; ++i_z) {
                int z_left = i_z-1;
                int z_right = i_z+1;
                if (z_left == -1) {
                    z_left = shape[1]-1;
                }
                if (z_right == shape[1]) {
                    z_right = 0;
                }
                double dphi_dr = (phi(r_right,i_theta,i_z) - 
                            phi(r_left,i_theta,i_z))/r_step;
                double dphi_dtheta = (phi(i_r,theta_right,i_z) - 
                            phi(i_r,theta_left,i_z))/theta_step;
                double dphi_dz = (phi(i_r,i_theta,z_right) - 
                            phi(i_r,i_theta,z_left))/z_step;
                Ex(i_r,i_theta,i_z) = cos(theta)*dphi_dr - sin(theta)*dphi_dtheta/r;
                Ey(i_r,i_theta,i_z) = sin(theta)*dphi_dr + cos(theta)*dphi_dtheta/r;
                Ez(i_r,i_theta,i_z) = dphi_dz;
            }
        }
    }
}

/*void
apply_E_n_kick_cylindrical(const Field_domain &fdomain,
                           Array_3d<double> &E, int n, double tau,
                           Macro_bunch_store &mbs)
{
    //jfa stub
}*/

inline
double
interpolate_3d(double x1, double x2, double x3,
               const Field_domain &fdomain,
               const Array_3d<double> &points)
{
    // Interpolate between grid points. There is no unique scheme to do this
    // in 3D, so we choose to use trilinear interpolation.
    std::vector<int> c(3); // c for corner
    std::vector<double> f(3); // f for fractional difference
    fdomain.get_leftmost_indices_offsets(x1,x2,x3,c,f);
    double val = ((1.0 - f[0]) * (1.0 - f[1]) * (1.0 - f[2]) * points(c[0],c[1],c[2]) +
            f[0] * (1.0 - f[1]) * (1.0 - f[2]) * points(c[0] + 1, c[1], c[2]) +
            (1.0 - f[0]) * f[1] * (1.0 - f[2]) * points(c[0], c[1] + 1, c[2]) +
            (1.0 - f[0]) * (1.0 - f[1]) * f[2] * points(c[0], c[1], c[2] + 1) +
            f[0] * f[1] * (1.0 - f[2]) * points(c[0] + 1, c[1] + 1, c[2]) +
            f[0] * (1.0 - f[1]) * f[2] * points(c[0] + 1, c[1], c[2] + 1) +
            (1.0 - f[0]) * f[1] * f[2] * points(c[0], c[1] + 1, c[2] + 1) +
            f[0] * f[1] * f[2] * points(c[0] + 1, c[1] + 1, c[2] + 1));
    return val;
}

void
full_kick_cylindrical(const Field_domain &fdomain,
                      Array_3d<double> &phi, double tau, 
                      Macro_bunch_store &mbs, Array_2d<double> &coords)
{
    std::vector<int> shape = fdomain.get_grid_shape();
    Array_3d<double> Ex(shape[0],shape[1],shape[2]);
    Array_3d<double> Ey(shape[0],shape[1],shape[2]);
    Array_3d<double> Ez(shape[0],shape[1],shape[2]);
    calculate_E_cylindrical(fdomain,phi,Ex,Ey,Ez);
    
    // jfa: I am taking this calculation of "factor" more-or-less
    // directly from Impact.  I should really redo it in terms that make
    // sense to me
    double gamma = -1 * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const  double c = 2.99792458e8;
    const  double pi = 4.0 * atan(1.0);

    double mass = mbs.mass * 1.0e9;
    double eps0 = 1.0 / (4 * pi * c * c * 1.0e-7); // using c^2 = 1/(eps0 * mu0)
    double Brho = gamma * beta * mass / c;
    double perveance0 = mbs.total_current / (2 * pi * eps0 * Brho * gamma * gamma*\
                beta * c * beta * c);
    double xk = mbs.units(0);

    double factor = pi * perveance0 * gamma * beta * beta / xk * 4.0 * pi;
    factor *= 1.0 / mbs.total_num;
    double zfactor = -tau* factor * beta * gamma * gamma;
    // (We think) this is for the Lorentz transformation of the transverse
    // E field.
    double xyfactor = -tau* factor * gamma;
    double jfa_total_kick_x = 0.0, jfa_total_kick_y = 0.0, jfa_total_kick_z = 0.0;
    for (int n = 0; n < mbs.local_num; ++n) {
        double r = coords(0,n);
        double theta = coords(1,n);
        double z = coords(2,n);
        //~ std::cout << "jfa: " << xyfactor << " " << interpolate_3d(r,theta,z,fdomain,Ex) << std::endl;
        mbs.local_particles(0,n) += xyfactor*
                interpolate_3d(r,theta,z,fdomain,Ex);
        jfa_total_kick_x += fabs(xyfactor*
                interpolate_3d(r,theta,z,fdomain,Ex));
        mbs.local_particles(2,n) += xyfactor*
                interpolate_3d(r,theta,z,fdomain,Ey);
        jfa_total_kick_y += fabs(xyfactor*
                interpolate_3d(r,theta,z,fdomain,Ey));
        mbs.local_particles(4,n) += zfactor*
                interpolate_3d(r,theta,z,fdomain,Ez);
        jfa_total_kick_z += fabs(zfactor*
                interpolate_3d(r,theta,z,fdomain,Ez));
    }
    std::cout << "jfa average kick x: " << jfa_total_kick_x/mbs.local_num << std::endl;
    std::cout << "jfa average kick y: " << jfa_total_kick_y/mbs.local_num << std::endl;
    std::cout << "jfa average kick z: " << jfa_total_kick_z/mbs.local_num << std::endl;
}
