#include "electric_field.h"
#include "communicate.h"
#include "mpi.h"
#include "mytimer.h"
#include <cmath>
#include "math_constants.h"
#include <fstream>
#include <stdio.h>

std::ofstream * fdebug=0;

void
init_fdebug()
{
    if (fdebug == 0) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        char buffer[100];
        sprintf(buffer,"cxxdebug-%02d",rank);
        fdebug = new std::ofstream(buffer);
    }
}
  
Real_scalar_field
calculate_E_n(Real_scalar_field &phi, int n)
{
    //~ *fdebug << "in calculate_E_n\n"; fdebug->flush();
    reset_timer();
    if ((n < 0) || (n > 2)) {
        std::stringstream message("");
        message << "calculate_E_n: invalid argument n=" << n
        << ". Argument be in range 0<=n<=2";
        throw std::invalid_argument(message.str());
    }
    Real_scalar_field E(phi.get_points().get_shape(), phi.get_physical_size(),
                        phi.get_physical_offset());
    E.get_points().zero_all();
    timer("E zero");
    Int3 shape(phi.get_points().get_shape());
    Double3 hi(phi.get_cell_size());
    double h(hi[n]), delta;
    Int3 point;
    Int3 offset_plus(0, 0, 0), offset_minus(0, 0, 0);
    offset_plus[n] = 1;
    offset_minus[n] = -1;
    double deriv;
    // calculate i range taking into account guard grids
    int i_lower, i_upper;
    i_lower = phi.get_points().get_dim0_lower();
    if (i_lower > 0) {
        i_lower += 1;
    }
    i_upper = phi.get_points().get_dim0_upper();
    if (i_upper < shape[0] - 1) {
        i_upper -= 1;
    }
    //~ *fdebug << "calculate_E_n loop is ilower iupper shape[1] shape[2] " 
        //~ << i_lower<<" "
        //~ << i_upper<<" "
        //~ << shape[1]<<" "
        //~ << shape[2]<<" "
        //~ << std::endl;
    //~ fdebug->flush();

    for (int i = i_lower; i < i_upper; ++i) {
        point[0] = i;
        for (int j = 0; j < shape[1]; ++j) {
            point[1] = j;
            for (int k = 0; k < shape[2]; ++k) {
                point[2] = k;
                Int3 p0(point), p1(point);
                if (point[n] == 0) {
                    p1.add(offset_plus);
                    delta = h;
                } else if (point[n] == shape[n] - 1) {
                    p0.add(offset_minus);
                    delta = h;
                } else {
                    p0.add(offset_minus);
                    p1.add(offset_plus);
                    delta = 2.0 * h;
                }
                deriv = (phi.get_points().get(p1) - phi.get_points().get(p0)) / delta;
                E.get_points().set(point, deriv);
            }
        }
    }
    timer("E calc");
    //~ *fdebug << "about to broadcast_E\n"; fdebug->flush();
    broadcast_E(E, i_lower, i_upper);
    //~ *fdebug << "broadcast_E done\n"; fdebug->flush();
    timer("E broadcast");
    return E;
}

void
apply_E_n_kick(Real_scalar_field &E, int n_axis, double tau,
               Macro_bunch_store &mbs)
{
    if ((n_axis < 0) || (n_axis > 2)) {
        std::stringstream message("");
        message << "apply_E_n_kick: invalid argument n_axis=" << n_axis
        << ". Argument be in range 0<=n_axis<=2";
        throw std::invalid_argument(message.str());
    }
    // jfa: I am taking this calculation of "factor" more-or-less
    // directly from Impact.  I should really redo it in terms that make
    // sense to me
    double gamma = -1 * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const  double c = 2.99792458e8;

    double mass = mbs.mass * 1.0e9;
    double eps0 = 1.0 / (4 * pi * c * c * 1.0e-7); // using c^2 = 1/(eps0 * mu0)
    double Brho = gamma * beta * mass / c;
    double perveance0 = mbs.total_current / (2 * pi * eps0 * Brho * gamma * gamma*\
                        beta * c * beta * c);
    double xk = mbs.units(0);

    double factor = pi * perveance0 * gamma * beta * beta / xk * 4.0 * pi;
    factor *= 1.0 / mbs.total_num;
    if (n_axis == 2) {
        factor *= beta * gamma * gamma;
    } else {
        // (We think) this is for the Lorentz transformation of the transverse
        // E field.
        factor *= gamma;
    }
    int index = 2 * n_axis + 1; // for axis n_axis = (0,1,2) corresponding to x,y,z,
    // in particle store indexing, px,py,pz = (1,3,5)
    double kick;
//     double jfa_total_kick = 0;
    for (int n = 0; n < mbs.local_num; ++n) {
        //~ if (n == 0) std::cout << "jfa: " << tau*factor << " " << E.get_val(Double3(mbs.local_particles(0, n),
            //~ mbs.local_particles(2, n),
                                //~ mbs.local_particles(4, n))) << std::endl;
        kick = tau * factor * E.get_val(Double3(mbs.local_particles(0, n),
                                                mbs.local_particles(2, n),
                                                mbs.local_particles(4, n)));
        mbs.local_particles(index, n) -= kick;
//         jfa_total_kick += fabs(kick);
    }
//     std::cout << "jfa: average kick for " << n_axis << ": " << jfa_total_kick/mbs.local_num << std::endl;
    timer("apply kick");
}

void
apply_phi_kick(Real_scalar_field &phi, int axis, double tau,
               Macro_bunch_store &mbs)
{
    if ((axis < 0) || (axis > 2)) {
        std::stringstream message("");
        message << "apply_E_n_kick: invalid argument axis=" << axis
        << ". Argument be in range 0<=axis<=2";
        throw std::invalid_argument(message.str());
    }
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
    if (axis == 2) {
        factor *= beta * gamma * gamma;
    } else {
        // (We think) this is for the Lorentz transformation of the transverse
        // E field.
        factor *= gamma;
    }
    int index = 2 * axis + 1; // for axis axis = (0,1,2) corresponding to x,y,z,
    // in particle store indexing, px,py,pz = (1,3,5)
    double kick;
    for (int n = 0; n < mbs.local_num; ++n) {
        kick = tau * factor * phi.get_deriv(Double3(mbs.local_particles(0, n),
                                            mbs.local_particles(2, n),
                                            mbs.local_particles(4, n)), axis);
        mbs.local_particles(index, n) -= kick;
    }
    timer("apply kick");
}

void
full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs)
{
    //~ init_fdebug();
    for (int axis = 0; axis < 3; ++axis) {
        //~ *fdebug << "about to Real_scalar_field\n"; fdebug->flush();
        Real_scalar_field E = calculate_E_n(phi, axis);
        //~ *fdebug << "about to apply kick " << axis << "\n"; fdebug->flush();
        apply_E_n_kick(E, axis, tau, mbs);
        //~ *fdebug << "full_kick complete\n"; fdebug->flush();
    }
}

void
transverse_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs)
{
    for (int axis = 0; axis < 2; ++axis) {
        Real_scalar_field E = calculate_E_n(phi, axis);
        apply_E_n_kick(E, axis, tau, mbs);
    }
}

void
just_phi_full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs)
{
    timer("misc");
    Real_scalar_field full_phi(phi.get_points().get_shape(), phi.get_physical_size(),
                               phi.get_physical_offset());
    allgather_phi(phi, full_phi);
    timer("E broadcast");
    for (int axis = 0; axis < 3; ++axis) {
        apply_phi_kick(full_phi, axis, tau, mbs);
        //~ apply_E_n_kick(E, axis, tau, mbs);
    }
}

void
rw_kick(Real_scalar_field &rho,
                Array_1d<double> &zdensity,
                Array_1d<double> &xmom, 
                Array_1d<double> &ymom,
                double tau, 
                Macro_bunch_store &mbs,
                double pipe_radiusx,
                double pipe_radiusy,
                double pipe_conduct,
                double zoffset)
{
    double gamma = -1 * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const double c = 2.99792458e8;
    double w = mbs.units(0)*c;

    double qe = 1.602176462e-19;
    double mass = mbs.mass * 1.0e9 * qe/(c*c); // convert mass in GeV to kg
    double eps0 = 1.0 / (4 * pi * c * c * 1.0e-7); // using c^2 = 1/(eps0 * mu0)

    double r_classical = mbs.charge*mbs.charge*qe*qe/(4*pi*eps0*mass*c*c);

    // tau is the step length in m
    double L = tau;
    
    // Number of particles in slice: Ntot_real*N_macro(slice)/Ntot_macro
    // N_macro(slice) is simply zdensity
    // Ntot_macro is mbs.total_num
    double Qtot = 2*pi*gamma*mbs.total_current/w;
    double Ntot_real = Qtot/(mbs.charge*qe);
    // N = N_factor * N_macro(slice) = N_factor * zdensity(slice)
    double N_factor = Ntot_real/mbs.total_num;
    
    // formula from paper is for delta pxy/p. We need the change
    // in trans mom coord, delta pxy/(mc)
    double dpop_to_trans_coord_factor = gamma*beta;
    
    double wake_factor_x = r_classical*2/
    		(beta*gamma*pi*pipe_radiusx*pipe_radiusx*pipe_radiusx)*
    		sqrt(4*pi*eps0*c/pipe_conduct)*L*
    		N_factor*dpop_to_trans_coord_factor;
    double wake_factor_y = r_classical*2/
    		(beta*gamma*pi*pipe_radiusy*pipe_radiusy*pipe_radiusy)*
    		sqrt(4*pi*eps0*c/pipe_conduct)*L*
    		N_factor*dpop_to_trans_coord_factor;
    
    int num_slices = zdensity.get_shape()[0];
    double left_z = rho.get_left()[2];
    double cell_size_z = rho.get_cell_size()[2];
    for (int n = 0; n < mbs.local_num; ++n) {
        int first_ahead_slice;
        if (zoffset == 0.0) {
            first_ahead_slice = static_cast<int>
                            (floor((mbs.local_particles(4,n) - left_z) 
                            / cell_size_z))+1;
            if (first_ahead_slice < 0) {
            	first_ahead_slice = num_slices;
            }
        } else {
            first_ahead_slice = 0;
        }
        double xkick, ykick;
        xkick = 0.0;
        ykick = 0.0;
        for (int ahead_slice = first_ahead_slice; 
                ahead_slice < num_slices;
                ++ahead_slice) {
            double zdistance_beamframe = (ahead_slice+0.5)*cell_size_z+left_z -
                mbs.local_particles(4,n)+zoffset*gamma;
            double zdistance = zdistance_beamframe/gamma;
            if (zdistance>0.0) {
                xkick += wake_factor_x * zdensity(ahead_slice) *
                    xmom(ahead_slice)/sqrt(zdistance);
                ykick += wake_factor_y * zdensity(ahead_slice) * 
                    ymom(ahead_slice)/sqrt(zdistance);
            } else {
                std::cerr << "warning: rw_kick encountered a nonsensical longitudinal distance\n";
            }

        }
        mbs.local_particles(1,n) += xkick;
        mbs.local_particles(3,n) += ykick;
    }
}

