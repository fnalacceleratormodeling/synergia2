#include "electric_field.h"
#include "communicate.h"
#include "mpi.h"
#include "mytimer.h"
#include "math_constants.h"

Real_scalar_field
calculate_E_n(Real_scalar_field &phi, int n)
{
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
    broadcast_E(E, i_lower, i_upper);
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
    for (int n = 0; n < mbs.local_num; ++n) {
        //~ if (n == 0) std::cout << "jfa: " << tau*factor << " " << E.get_val(Double3(mbs.local_particles(0, n),
            //~ mbs.local_particles(2, n),
                                //~ mbs.local_particles(4, n))) << std::endl;
        kick = tau * factor * E.get_val(Double3(mbs.local_particles(0, n),
                                                mbs.local_particles(2, n),
                                                mbs.local_particles(4, n)));
        mbs.local_particles(index, n) -= kick;
    }
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
    for (int axis = 0; axis < 3; ++axis) {
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
    // jfa: I am taking this calculation of "factor" more-or-less
    // directly from Impact.  I should really redo it in terms that make
    // sense to me
    double gamma = -1 * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const  double c = 2.99792458e8;

    double mass = mbs.mass * 1.0e9;
    double eps0 = 1.0 / (4 * pi * c * c * 1.0e-7); // using c^2 = 1/(eps0 * mu0)
    double qe = 1.602176462e-19;
    double emradius = qe/(4*pi*eps0*mass);

    // jfa: from Eric's Impedance.f90
            //~ wake_factor = 2.0/(pi * pipe_radius**3) *  &
             //~ sqrt(speed_of_light/(pipe_conduct/fourpiepsilon0)) * &
             //~ accel_length

    double accel_length = tau;
    double wake_factorx = 2.0/(pi*pipe_radiusx*pipe_radiusx*pipe_radiusx)*
        sqrt(c/(pipe_conduct/(4*pi*eps0)))*accel_length;
    double wake_factory= 2.0/(pi*pipe_radiusy*pipe_radiusy*pipe_radiusy)*
        sqrt(c/(pipe_conduct/(4*pi*eps0)))*accel_length;
    
    double charge_factor = mbs.total_current/mbs.total_num;
    
    int num_slices = zdensity.get_shape()[0];
    double left_z = rho.get_left()[2];
    double cell_size_z = rho.get_cell_size()[2];
    for (int n = 0; n < mbs.local_num; ++n) {
        int first_ahead_slice;
        if (zoffset == 0.0) {
            first_ahead_slice = static_cast<int>
                            (floor((mbs.local_particles(4,n) - left_z) 
                            / cell_size_z))+1;
        } else {
            first_ahead_slice = 0;
        }
        double xkick, ykick;
        xkick = 0.0;
        ykick = 0.0;
        for (int ahead_slice = first_ahead_slice; 
                ahead_slice < num_slices;
                ++ahead_slice) {
            double zdistance = (ahead_slice+0.5)*cell_size_z+left_z -
                mbs.local_particles(4,n)+zoffset;
            if (zdistance>0.0) {
                xkick += wake_factorx * charge_factor* zdensity(ahead_slice) *
                    emradius/(beta*gamma)*
                    xmom(ahead_slice)/sqrt(zdistance);
                ykick += wake_factory * charge_factor * zdensity(ahead_slice) * 
                    emradius/(beta*gamma)*
                    ymom(ahead_slice)/sqrt(zdistance);
            } else {
                //~ std::cout << "jfa: zdistance = " << zdistance << std::endl;
                //~ std::cout << "jfa: zoffset = " << zoffset << std::endl;
                //~ std::cout << "jfa: cell_size_z = " << cell_size_z << std::endl;
                //~ std::cout << "jfa: left_z = " << left_z << std::endl;
                //~ std::cout << "jfa: ahead_slice = " << ahead_slice << std::endl;
                //~ exit(1);
            }

        }
        mbs.local_particles(1,n) += xkick;
        mbs.local_particles(3,n) += ykick;
    }
}

