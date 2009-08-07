#include "electric_field.h"
#include "communicate.h"
#include "mpi.h"
#include "mytimer.h"
#include <cmath>
#include "math_constants.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "basic_toolkit/PhysicsConstants.h"



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
calculate_E_n(Real_scalar_field &phi, int n, Fftw_helper &fftwh, bool z_periodic)
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

     int i_lower, i_upper;   
     i_lower =fftwh.lower(); 
     i_upper =std::min(fftwh.upper(),shape[0]);


// AM: the field on  the grid edges is set to zero for directions with no periodic boundary, 
//in acordance  to charge deposit, thus  we loop from 1 to shape-2 
     int i_lower1, i_upper1; 
     i_lower1=std::max(i_lower,1);
     i_upper1=std::min(i_upper,shape[0]-1);

    for (int i = i_lower1; i < i_upper1; ++i) {
        point[0] = i;
        for (int j = 1; j < shape[1]-1; ++j) { // the field on  the grid edges is set to zero
            point[1] = j;
            for (int k = 0; k < shape[2]; ++k) { // z direction may be periodic
                point[2] = k;
                Int3 p0(point), p1(point);
                if (point[n] == 0) {
                    p1.add(offset_plus);
                    delta = h;
                    deriv=0.;
                } else if (point[n] == shape[n] - 1) {
                    p0.add(offset_minus);
                    delta = h;
                    deriv=0.;
                } else {
                    p0.add(offset_minus);
                    p1.add(offset_plus);

                    delta = 2.0 * h;
                }
                deriv = (phi.get_points().get(p1) - phi.get_points().get(p0))/ delta; 
                E.get_points().set(point, deriv);
	        
            }
        }
    }
   
    
      
    for (int i = i_lower1; i < i_upper1; ++i) {
      for (int j = 0; j < shape[1]; ++j) { 
         Int3 left(i, j, 0), right(i, j, shape[2] - 1);
         double average;
         if ((z_periodic))  {
                   average=0.5*(E.get_points().get(left)+E.get_points().get(right));
          } else {average=0.;} // the field on  the grid edges is set to zero
           E.get_points().set(left, average);
 	   E.get_points().set(right, average);
      } 
     }
      
/*
    for (int j = 0; j < shape[1]; ++j) {
      for (int k = 0; k < shape[2]; ++k) { 
	Int3 left(0, j, k), right(shape[0]-1, j, k);	
 	E.get_points().set(left, 0.);
 	E.get_points().set(right, 0.);
      }
    }

    for (int i = 0; i < shape[0]; ++i) {
      for (int k = 0; k < shape[2]; ++k) { 
	Int3 left(i, 0, k), right(i, shape[1]-1, k);	
 	E.get_points().set(left, 0.);
 	E.get_points().set(right, 0.);
      }
    }*/


  
	
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


 
    double gamma = -1. * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const  double c = PH_MKS_c;
    double mass = mbs.mass * 1.0e9;
    double eps0 = PH_MKS_eps0; 

/*    consider the electric fields in the rest frame, E'x, E'y, E'z,
	and Ex,Ey,Ez in the lab frame
      the force on a charge q is
	Fx=q * E'x/gamma =q * Ex/gamma^2
	Fy=q * E'y/gamma=q * Ey/gamma^2
	Fz=q * E'z = q *Ez	
	
	the kick in the 3rd direction is a kick of p_t, not p_z	
	p_t=-U ==> dp_t/dt=-betaz * dpz/dt-betax *dpx/dt-betay*dpy/dt =~ -beta*dpz/dt
*/

    double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
    double factor =PH_CNV_brho_to_p/eps0; // charge of the particle is PH_CNV_brho_to_p =p/Brho in [GeV/c]/[Tesla*m]
    factor *= length* mbs.total_current /(beta * c); // total charge=linear charge density*length
    factor *= 1.0/(beta * c); // the  arc length tau=beta*c* (Delta t), so (Delta t)= tau/(beta*c)
    factor *= -1.0/gamma;    // relativistic factor,...-minus from the definition of E= + grad phi ???
    factor *=mbs.units(1); // the kikcing force should be muliplied  by the unit of p, this is a factor of 1/mass
    if (n_axis == 2) {factor *=-beta*gamma;} // -dp_t=-beta dp_z; E      
    int index = 2 * n_axis + 1; // for axis n_axis = (0,1,2) corresponding to x,y,z,
    // in particle store indexing, px,py,pz = (1,3,5)
    double kick;
/*   difference with the factor (i.e. xycon,tcon) in impact :
 	   factor_fish*n_part/(4*pi)= - factor_impact
 
      in impact, -1/(4*pi) *n_part is included in the definition of the electric field....
*/
	
    for (int n = 0; n < mbs.local_num; ++n) {
       
        kick = tau * factor * E.get_val(Double3(mbs.local_particles(0, n),
                                                mbs.local_particles(2, n),
                                                mbs.local_particles(4, n)));

   
        mbs.local_particles(index, n) += kick;
	
    }

   	


    timer("apply kick");
	
}



void apply_Efield_kick(const std::vector<Real_scalar_field> &E, double tau,
               Macro_bunch_store &mbs)
{
 
    double gamma = -1. * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const  double c = PH_MKS_c;
    double mass = mbs.mass * 1.0e9;
    double eps0 = PH_MKS_eps0; 

    double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
    double factor =PH_CNV_brho_to_p/eps0; // charge of the particle is PH_CNV_brho_to_p =p/Brho
    factor *= length* mbs.total_current /(beta * c); // total charge=linear charge density*length
    factor *= 1.0/(beta * c); // the  arc length tau=beta*c* (Delta t), so (Delta t)= tau/(beta*c)
    factor *= -1.0/gamma;    // relativistic factor,...-minus from the definition of E= + grad phi ???
    factor *=mbs.units(1); // the kikcing force should be muliplied  by the unit of p, this is a factor of 1/mass



/*    consider the electric fields in the rest frame, E'x, E'y, E'z,
	and Ex,Ey,Ez in the lab frame
      the force on a charge q is
	Fx=q * E'x/gamma =q * Ex/gamma^2
	Fy=q * E'y/gamma=q * Ey/gamma^2
	Fz=q * E'z = q *Ez	
	
	the kick in the 3rd direction is a kick of p_t, not p_z	
	p_t=-U ==> dp_t/dt=-betaz * dpz/dt-betax *dpx/dt-betay*dpy/dt =~ -beta*dpz/dt
*/

/*   difference with the factor (i.e. xycon,tcon) in impact :
 	   factor_fish*n_part/(4*pi)= - factor_impact
 
 in impact, -1/(4*pi) *n_part is included in the definition of the electric field....
*/
    double kick;
    double factor1;
    int index;
    for (int axis = 0; axis < 3; ++axis) {
          index = 2 * axis + 1;
          axis == 2 ? factor1 =-factor*beta*gamma:factor1=factor;
          for (int n = 0; n < mbs.local_num; ++n) {                     
                             kick = tau * factor1* E.at(axis).get_val(Double3(mbs.local_particles(0, n),
                                               mbs.local_particles(2, n),
                                               mbs.local_particles(4, n)));
              mbs.local_particles(index, n) += kick;
           }
    }
  return;
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


     double gamma = -1. * mbs.ref_particle(5);
     double beta = sqrt(gamma * gamma - 1.0) / gamma;
     const  double c = PH_MKS_c;
     double mass = mbs.mass * 1.0e9;
     double eps0 = PH_MKS_eps0; 

/*    consider the electric fields in the rest frame, E'x, E'y, E'z,
	and Ex,Ey,Ez in the lab frame
      the force on a charge q is
	Fx=q * E'x/gamma =q * Ex/gamma^2
	Fy=q * E'y/gamma=q * Ey/gamma^2
	Fz=q * E'z = q *Ez	
	
	the kick in the 3rd direction is a kick of p_t, not p_z	
	p_t=-U ==> dp_t/dt=-betaz * dpz/dt-betax *dpx/dt-betay*dpy/dt =~ -beta*dpz/dt
*/

    double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
    double factor =PH_CNV_brho_to_p/eps0; // charge of the particle is PH_CNV_brho_to_p =p/Brho
    factor *= length* mbs.total_current /(beta * c); // total charge=linear charge density*length
    factor *= 1.0/(beta * c); // the  arc length tau=beta*c* (Delta t), so (Delta t)= tau/(beta*c)
    factor *= -1.0/gamma;    // relativistic factor...,-minus from the definition of E= + grad phi ???
    factor *=mbs.units(1); // the kikcing force should be muliplied  by the unit of p, this is a factor of 1/mass
    if (axis == 2) {factor *=-beta*gamma;} // -dp_t=-beta dp_z; E     


/*   difference with the factor (i.e. xycon,tcon) in impact :
 	   factor_fish*n_part/(4*pi)= - factor_impact
 
 in impact, -1/(4*pi) *n_part is included in the definition of the electric field....
*/
  
 
    int index = 2 * axis + 1; // for axis n_axis = (0,1,2) corresponding to x,y,z,
    // in particle store indexing, px,py,pz = (1,3,5)
    double kick;   
    for (int n = 0; n < mbs.local_num; ++n) {
        kick = tau * factor * phi.get_deriv(Double3(mbs.local_particles(0, n),
                                            mbs.local_particles(2, n),
                                            mbs.local_particles(4, n)), axis);
        mbs.local_particles(index, n) += kick;	

    }
    timer("apply kick");
      
}

void
full_kick_version(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs,Fftw_helper &fftwh, bool z_periodic)
{
    //~ init_fdebug();
   std::vector<Real_scalar_field> E;	
    for (int axis = 0; axis < 3; ++axis) {
        //~ *fdebug << "about to Real_scalar_field\n"; fdebug->flush();
	Real_scalar_field En=calculate_E_n(phi, axis, fftwh, z_periodic);
        E.push_back(En);
           }
        //~ *fdebug << "about to apply kick " << axis << "\n"; fdebug->flush();
        apply_Efield_kick(E, tau, mbs);
        //~ *fdebug << "full_kick complete\n"; fdebug->flush();
   
}


void
full_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs, Fftw_helper &fftwh, bool z_periodic )
{
    //~ init_fdebug();
    for (int axis = 0; axis < 3; ++axis) {
        //~ *fdebug << "about to Real_scalar_field\n"; fdebug->flush();
        Real_scalar_field E = calculate_E_n(phi, axis, fftwh, z_periodic);
        //~ *fdebug << "about to apply kick " << axis << "\n"; fdebug->flush();
        apply_E_n_kick(E, axis, tau, mbs);
        //~ *fdebug << "full_kick complete\n"; fdebug->flush();
     //   apply_phi_kick(phi, axis,  tau,  mbs);
    }
}


void
transverse_kick(Real_scalar_field &phi, double tau, Macro_bunch_store &mbs, Fftw_helper &fftwh, bool z_periodic)
{
    for (int axis = 0; axis < 2; ++axis) {
        Real_scalar_field E = calculate_E_n(phi, axis, fftwh, z_periodic);
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
rw_kick_eric(double left_z, double size_z,
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

    // XXXXXXXXXXXXXXX  dbg is set here
    int dbg = 0;
    
    if (dbg) std::cout << "rwkick called, w: " << w <<
      " total current: " << mbs.total_current <<
      " L: " << L << std::endl;

    if (dbg) std::cout << "left: " << left_z << " size: " << size_z << std::endl;

    // Number of particles in slice: Ntot_real*N_macro(slice)/Ntot_macro
    // N_macro(slice) is simply zdensity
    // Ntot_macro is mbs.total_num
    double Qtot = 2*pi*mbs.total_current/w;
    double Ntot_real = Qtot/(mbs.charge*qe);
    // N = N_factor * N_macro(slice) = N_factor * zdensity(slice)
    double N_factor = Ntot_real/mbs.total_num;
    

    // formula from paper is for delta pxy/p. We need the change
    // in trans mom coord, delta pxy/(mc)
    double dpop_to_trans_coord_factor = gamma*beta;
    
    if (dbg) std::cout << "is_fixed_z: " << mbs.is_fixedz << " mbs.total_current: " <<
      mbs.total_current << std::endl;
    if (dbg) std::cout << "Qtot: " << Qtot << "Ntot_real: " << Ntot_real <<
      "total_num: " << mbs.total_num << "beta*gamma: " <<
	       dpop_to_trans_coord_factor << " N_factor: " << N_factor << std::endl;

    if (dbg) std::cout << "pipe_conduct: " << pipe_conduct << std::endl;

    double wake_factor_x = r_classical*2/
    		(beta*gamma*pi*pipe_radiusx*pipe_radiusx*pipe_radiusx)*
    		sqrt(4*pi*eps0*c/pipe_conduct)*L*
    		N_factor * dpop_to_trans_coord_factor;
    double wake_factor_y = r_classical*2/
    		(beta*gamma*pi*pipe_radiusy*pipe_radiusy*pipe_radiusy)*
    		sqrt(4*pi*eps0*c/pipe_conduct)*L*
                N_factor * dpop_to_trans_coord_factor;
    
    if (dbg) std::cout << "wake_factor_x: " << wake_factor_x << " wake_factor_y: " <<
      wake_factor_y << std::endl;
    int num_slices = zdensity.get_shape()[0];
    double cell_size_z = size_z/num_slices;

    if (dbg) std::cout << "num_slices: " << num_slices << "cell_size_z: " << cell_size_z << std::endl;
    if (dbg) {
      double sumdens = 0.0;
      for (int k=0; k<num_slices; ++k)
	sumdens += zdensity(k);
      if (dbg) std::cout << "Total density sum: " << sumdens << std::endl;
    }

    for (int n = 0; n < mbs.local_num; ++n) {
        int first_ahead_slice;
        if (zoffset == 0.0) {
            first_ahead_slice = static_cast<int>
                            (floor((mbs.local_particles(4,n) - left_z) 
                            / cell_size_z))+1;
            if (first_ahead_slice < 0) {
            	first_ahead_slice = num_slices;
            }
	    if (((n % 1000) == 0) && (n < 100001) && (n > 0)) {
	      if (dbg) std::cout << "particle: " << n << " first_ahead_slice: " <<
		first_ahead_slice << std::endl;
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
	    if (((n % 1000) == 0) && (n < 100001) && (n > 0)) {
	      if (dbg) std::cout << "ahead_slice: " << ahead_slice <<
		" zdistance: " << zdistance <<
			 " zdensity: " << zdensity(ahead_slice) <<
			 " xmom, ymom: " <<
			 xmom(ahead_slice) << " " << ymom(ahead_slice) << std::endl;

	      if (zdistance>0.0) {
                if (dbg) std::cout << "x_kick: " << wake_factor_x * zdensity(ahead_slice) *
		  xmom(ahead_slice)/sqrt(zdistance) << " y_kick: " <<
                wake_factor_y * zdensity(ahead_slice) * 
			   ymom(ahead_slice)/sqrt(zdistance) << std::endl << std::endl;
	      } else {
                std::cerr << "warning: rw_kick encountered a nonsensical longitudinal distance\n";
	      }
	    }
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
    if (dbg) std::cout << std::endl;
}




void
rw_kick(        double size_z,
                Array_1d<int> &bin_partition,
                Array_1d<double> &zdensity,
                Array_1d<double> &xmom, 
                Array_1d<double> &ymom,
                double tau, 
                Macro_bunch_store &mbs,
                double w_f,
                double cutoff_small_z, Array_1d<double> &wake_coeff, 
                double quad_wake_sum, bool bool_quad_wake)

{


//  general  form of kikcks
//       xkick=wake_factor*(ax_dipole*dipole_x(bin)+
//                     bx_dipole*dipole_y(bin)+
//                    (cx_quad+a_quad*mbs.local_particles(0,n)+b_quad*mbs.local_particles(2,n))*quad(bin));
// 
//       ykick=wake_factor*(ay_dipole*dipole_y(bin)+
//                    by_dipole*dipole_x(bin)+
//                    (cy_quad-a_quad*mbs.local_particles(2,n)+b_quad*mbs.local_particles(0,n))*quad(bin));

// some parameters are zero due to symmetry, see S. Heifets, SLAC/AP110, January 1998......
//input parameter (ax_dipole, ay_dipole, bx_dipole, by_dipole,a_quad, b_quad, cx_quad, cy_quad)


 
    double ax_dipole= wake_coeff(0);
    double ay_dipole= wake_coeff(1);
    double bx_dipole= wake_coeff(2);
    double by_dipole= wake_coeff(3);
    double a_quad   = wake_coeff(4);
    double b_quad   = wake_coeff(5);
    double cx_quad  = wake_coeff(6);
    double cy_quad  = wake_coeff(7);






    double gamma = -1 * mbs.ref_particle(5);
    double beta = sqrt(gamma * gamma - 1.0) / gamma;
    const double c = PH_MKS_c;
    double qe = PH_MKS_e; // 1.602176462e-19;
    double eps0 = PH_MKS_eps0;// 1.0 / (4 * pi * c * c * 1.0e-7); // using c^2 = 1/(eps0 * mu0)
    double r_classical = PH_MKS_rp;//mbs.charge*mbs.charge*qe*qe/(4*pi*eps0*mass*c*c);
    // tau is the step length in m
    double L = tau;

    
    // Number of particles in slice: Ntot_real*N_macro(slice)/Ntot_macro
    // N_macro(slice) is simply zdensity
    // Ntot_macro is mbs.total_num
   
    double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
    double Qtot = length* mbs.total_current /(beta * c);
    double Ntot_real = Qtot/(mbs.charge*qe);  
     //std::cout<<" Ntot_real ="<<Ntot_real<<std::endl;
    // N = N_factor * N_macro(slice) = N_factor * zdensity(slice)
    double N_factor = Ntot_real/mbs.total_num;
    
  

//    w_f= paking_fraction*r_classical*2.0/
//     		(beta*gamma*pi*pipe_radius*pipe_radius*pipe_radius)*
//     		sqrt(4*pi*eps0*c/pipe_conduct)*gamma*beta; the multiplication with gamma*beta of
//             formula from paper is to account that our coordinate is delta pxy/(mc) not delta pxy/p

    double  wake_factor=w_f*N_factor*L;

   
    int num_slices = zdensity.get_shape()[0];
    Array_1d<double> dipole_x(num_slices);
    Array_1d<double> dipole_y(num_slices);
    Array_1d<double> quad(num_slices);
    double cell_size_z = size_z/num_slices;
    int cut_scaled=static_cast<int>(floor(cutoff_small_z*gamma/cell_size_z));	
    //std::cout<<" cutoff small="<< cutoff_small_z*gamma<<" cell size="<< cell_size_z<<"  cut_scaled="<<cut_scaled<<std::endl;

    get_wake_factors(num_slices, cut_scaled, zdensity, xmom, ymom, dipole_x, dipole_y, quad);


    wake_factor *= sqrt(gamma/cell_size_z); // the distance in lab frame is the distance 
                                          //  in the beam frame divided to gamma

    // contributions from previous turns, propto sum_n W(n*orbit_lenght), are considered by adding quad_wake_sum
    if (bool_quad_wake) { 
          double quad_wake_sum_scaled=mbs.total_num*quad_wake_sum*sqrt(cell_size_z/gamma); //rescaled to cancel
          quad.add(quad_wake_sum_scaled);                     //  the factor considered in wake_factor above
    }



//  applying kikcs	
   for (int n = 0; n < mbs.local_num; ++n) {
       double xkick=0., ykick=0.;
       int bin=bin_partition(n);  // bin_partition(n) is the bin where you find the particle n 
        /*   if ((bin>=num_slices) || (bin<0))  { std::cout<<"bin="<<bin<<"num_slices="<<num_slices<<std::endl;
  		                         throw
                                         std::runtime_error("something wrong with bining");}*/	


       xkick=ax_dipole*dipole_x(bin);
       ykick=ay_dipole*dipole_y(bin);
    //xkick=ax_dipole*dipole_x(bin)+bx_dipole*dipole_y(bin);
    //ykick=ay_dipole*dipole_y(bin)+by_dipole*dipole_x(bin);

       if (bool_quad_wake) {
                   xkick += a_quad*quad(bin)*mbs.local_particles(0,n);
                   ykick +=-a_quad*quad(bin)*mbs.local_particles(2,n);
      
              // xkick += (cx_quad+a_quad*mbs.local_particles(0,n)+b_quad*mbs.local_particles(2,n))*quad(bin);  //this is a more general form,
                                                                                                                // comment before and uncomment here
             //  ykick += (cy_quad-a_quad*mbs.local_particles(2,n)+b_quad*mbs.local_particles(0,n))*quad(bin);

       }	

      // mbs.local_particles(5,n) += -wake_factor*(xkick*mbs.local_particles(1,n)+ykick*mbs.local_particles(3,n))/gamma; //this is a higher order approximation
                                                                                                                         // probably not necessary...

       mbs.local_particles(1,n) += wake_factor*xkick;	
       mbs.local_particles(3,n) += wake_factor*ykick;
   }   
   
}

void get_wake_factors(int num_slices, int icut, Array_1d<double> &zdensity, 
Array_1d<double> &xmom, Array_1d<double> &ymom, Array_1d<double> &dipole_x, 
Array_1d<double> &dipole_y, Array_1d<double> & quad)
{

     dipole_x.set_all(0.0);
     dipole_y.set_all(0.0);
         quad.set_all(0.0);

     for (int i = 0; i < num_slices; ++i){      
       for (int j = i+1+icut; j < num_slices; ++j){  // icut is introduced to avoid the interaction at very small distance 
 	  dipole_x(i) += zdensity(j)*xmom(j)/sqrt(double(j-i)); // for resonable cell size (lgridnum<500) it is usually zero
          dipole_y(i) += zdensity(j)*ymom(j)/sqrt(double(j-i));
          quad(i) += zdensity(j)/sqrt(double(j-i)); 
       } 
 
   }


}







