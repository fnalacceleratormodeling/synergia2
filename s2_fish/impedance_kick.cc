#include "impedance_kick.h"
#include <cmath>
#include "math_constants.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "basic_toolkit/PhysicsConstants.h"




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
    //Qtot=mbs.bunch_np*qe*mbs.charge;// total charge
    double Ntot_real = Qtot/(mbs.charge*qe);  
     //std::cout<<" Ntot_real ="<<Ntot_real<<std::endl;
    // N = N_factor * N_macro(slice) = N_factor * zdensity(slice)
    double N_factor = Ntot_real/mbs.total_num;
    //double N_factor = mbs.bunch_np/mbs.total_num;

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
   // int cut_scaled=max(static_cast<int>(floor(cutoff_small_z*gamma/cell_size_z)),num_slices+1);	
   // std::cout<<" cutoff small="<< cutoff_small_z*gamma<<" cell size="<< cell_size_z<<"  cut_scaled="<<cut_scaled<<std::endl;
   // std::cout<<" gamma="<<gamma<<std::endl;
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


void rw_kick_transverse( double size_z,
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
    //Qtot=mbs.bunch_np*qe*mbs.charge;// total charge
    double Ntot_real = Qtot/(mbs.charge*qe);  
     //std::cout<<" Ntot_real ="<<Ntot_real<<std::endl;
    // N = N_factor * N_macro(slice) = N_factor * zdensity(slice)
    double N_factor = Ntot_real/mbs.total_num;
    //double N_factor = mbs.bunch_np/mbs.total_num;

    //    w_f= paking_fraction*r_classical*2.0/
//          (beta*gamma*pi*pipe_radius*pipe_radius*pipe_radius)*
//          sqrt(4*pi*eps0*c/pipe_conduct)*gamma*beta; the multiplication with gamma*beta of
//             formula from paper is to account that our coordinate is delta pxy/(mc) not delta pxy/p

    double  wake_factor=w_f*N_factor*L;

   
    int num_slices = zdensity.get_shape()[0];
    Array_1d<double> dipole_x(num_slices);
    Array_1d<double> dipole_y(num_slices);
    Array_1d<double> quad(num_slices);
    double cell_size_z = size_z/num_slices;
    int cut_scaled=static_cast<int>(floor(cutoff_small_z*gamma/cell_size_z));
   // int cut_scaled=max(static_cast<int>(floor(cutoff_small_z*gamma/cell_size_z)),num_slices+1);   
   // std::cout<<" cutoff small="<< cutoff_small_z*gamma<<" cell size="<< cell_size_z<<"  cut_scaled="<<cut_scaled<<std::endl;
   // std::cout<<" gamma="<<gamma<<std::endl;
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





