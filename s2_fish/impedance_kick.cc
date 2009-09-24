#include "impedance_kick.h"
#include <cmath>
#include "mpi.h"
#include "math_constants.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "basic_toolkit/PhysicsConstants.h"




void
rw_kick(       Array_1d<double> &dparameters,
                Array_1d<int> &bin_partition,
                Array_1d<double> &zdensity,
                Array_1d<double> &xmom, 
                Array_1d<double> &ymom,
                Macro_bunch_store &mbs,
                Array_1d<double> &wake_coeff, 
                bool bool_quad_wake, int bunch_i,
                Array_3d<double> &stored_means,
                Array_2d<int>  &stored_buckets,
                Array_2d<double>  &stored_bunchnp
                
        )
                

{
    
    double size_z= dparameters(0);
    double tau= dparameters(1);
    double w_f= dparameters(2);
    double cutoff_small_z=  dparameters(3);
    double quad_wake_sum=dparameters(4);
    double bunch_sp= dparameters(5);
    double line_length=dparameters(6) ; 
    
 
//  general  form of kikcks
//       xkick=wake_factor*(ax_dipole*dipole_x(bin)+
//                     bx_dipole*dipole_y(bin)+
//                    (cx_quad+a_quad*mbs.local_particles(0,n)+b_quad*mbs.local_particles(2,n))*quad(bin));
// 
//       ykick=wake_factor*(ay_dipole*dipole_y(bin)+
//                    by_dipole*dipole_x(bin)+
//                    (cy_quad-a_quad*mbs.local_particles(2,n)+b_quad*mbs.local_particles(0,n))*quad(bin));

// some parameters are zero due to symmetry, see S. Heifets, SLAC/AP110, January 1998......
//input parameter (ax_dipole, ay_dipole, bx_dipole, by_dipole,a_quad, b_quad, cx_quad, cy_quad, a_monopole=a0*b*b/4)


    double ax_dipole= wake_coeff(0);
    double ay_dipole= wake_coeff(1);
    double bx_dipole= wake_coeff(2);
    double by_dipole= wake_coeff(3);       // index 0 through 7, transverse impedance
    double a_quad   = wake_coeff(4);
    double b_quad   = wake_coeff(5);
    double cx_quad  = wake_coeff(6);
    double cy_quad  = wake_coeff(7);
    double a_monopole = wake_coeff(8);  // a0*pipe_radius^2/4 for longitudinal impedance






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
   
//     double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
//     double Qtot = length* mbs.total_current /(beta * c);
//     //Qtot=mbs.bunch_np*qe*mbs.charge;// total charge
//     double Ntot_real = Qtot/(mbs.charge*qe);  
//     std::cout<<" Ntot_real ="<<Ntot_real<<"  mbs.bunch_np="<<mbs.bunch_np<<std::endl;
    // N = N_factor * N_macro(slice) = N_factor * zdensity(slice)
//     double N_factor = Ntot_real/mbs.total_num;
    double N_factor = mbs.bunch_np/mbs.total_num;
    

//    w_f= paking_fraction*r_classical*2.0/
//     		(beta*gamma*pi*pipe_radius*pipe_radius*pipe_radius)*
//     		sqrt(4*pi*eps0*c/pipe_conduct)*gamma*beta; the multiplication with gamma*beta of
//             formula from paper is to account that our coordinate is delta pxy/(mc) not delta pxy/p

    double  wake_factor=w_f*N_factor*L;
    
   
    int num_slices = zdensity.get_shape()[0];
    Array_1d<double> dipole_x(num_slices);
    Array_1d<double> dipole_y(num_slices);
    Array_1d<double> quad(num_slices);
    Array_1d<double> l_monopole(num_slices);
    double cell_size_z = size_z/num_slices;
    double rescaling_f=gamma/cell_size_z;
    int cut_scaled=static_cast<int>(floor(cutoff_small_z*rescaling_f));
    double rbunch_sp=bunch_sp*rescaling_f; // rescaled bunch spacing
    
   // int cut_scaled=max(static_cast<int>(floor(cutoff_small_z*gamma/cell_size_z)),num_slices+1);	
   // std::cout<<" cutoff small="<< cutoff_small_z*gamma<<" cell size="<< cell_size_z<<"  cut_scaled="<<cut_scaled<<std::endl;
   // std::cout<<" gamma="<<gamma<<std::endl;
    get_wake_factors(num_slices, cut_scaled, zdensity, xmom, ymom, dipole_x, dipole_y, quad, l_monopole, 
                     stored_means,stored_buckets,stored_bunchnp, line_length, rbunch_sp,bunch_i,N_factor,rescaling_f);
 

    wake_factor *= sqrt(rescaling_f); // the distance in lab frame is the distance 
                                          //  in the beam frame divided to gamma
    
   
/*    the kick in the 3rd direction is a kick of p_t, not p_z 
            p_t=-U ==> dp_t/dt=-betaz * dpz/dt-betax *dpx/dt-betay*dpy/dt =~ -beta*dpz/dt */
     double  l_wake_factor=-beta*wake_factor*rescaling_f; // rescaling factor for z impedance is propto z^(3/2)
    
     // contributions from previous turns, propto sum_n W(n*orbit_lenght), are considered by adding quad_wake_sum
     // The assumption that bunch.total_num is not changing from previous tunrns is made!!!!
    if (bool_quad_wake) { 
// this is the contribution from turns > nstored turns
// if nstored turns  is large,  number of bunches can be taken  as a factor outside summation over bunches, 
//        i.e 1/(nL+bunch_sp) approx  1/(nL) , for n large
                   
        
        double quad_wake_sum_scaled=mbs.total_num*quad_wake_sum/sqrt(rescaling_f); //rescaled to cancel
                                                                               //  the factor considered in wake_factor above
        quad_wake_sum_scaled *= stored_buckets.get_shape()[0]; // =number of bunches....
          quad.add(quad_wake_sum_scaled); 
    }

    
//  applying kikcs	
   for (int n = 0; n < mbs.local_num; ++n) {
       double xkick=0., ykick=0., zkick=0.;
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
                                                                                                                // comment before and uncomment
                                                                                                                //here
             //  ykick += (cy_quad-a_quad*mbs.local_particles(2,n)+b_quad*mbs.local_particles(0,n))*quad(bin);

       }	
     
       zkick=a_monopole* l_monopole(bin);
     
                                                                                                                         // probably not necessary...
       
       mbs.local_particles(1,n) += wake_factor*xkick;	
       mbs.local_particles(3,n) += wake_factor*ykick;
       mbs.local_particles(5,n) += l_wake_factor*zkick;
        // mbs.local_particles(5,n) += -wake_factor*(xkick*mbs.local_particles(1,n)+ykick*mbs.local_particles(3,n))/gamma; //this is a higher order approximation
       
   }   
}

void get_wake_factors(int num_slices, int icut, Array_1d<double> &zdensity, 
Array_1d<double> &xmom, Array_1d<double> &ymom, Array_1d<double> &dipole_x, 
Array_1d<double> &dipole_y, Array_1d<double> & quad, Array_1d<double> & l_monopole, Array_3d<double> &stored_means,
Array_2d<int>  &stored_buckets, Array_2d<double>  &stored_bunchnp, double line_length, double rbunch_sp, 
int bunch_i, double N_factor, double rescaling_f)
{
     dipole_x.set_all(0.0);
     dipole_y.set_all(0.0);
     quad.set_all(0.0);
     l_monopole.set_all(0.0);
     
     
     int numbunches=stored_buckets.get_shape()[0];
     int registered_turns=stored_means.get_shape()[1];
     for (int i = 0; i < num_slices; ++i){   
       // in-bunch impedance     
       for (int j = i+1+icut; j < num_slices; ++j){  // icut is introduced to avoid the interaction at very small distance 
 	      dipole_x(i) += zdensity(j)*xmom(j)/sqrt(double(j-i)); // for resonable cell size (lgridnum<500) it is usually zero
          dipole_y(i) += zdensity(j)*ymom(j)/sqrt(double(j-i));
          quad(i) += zdensity(j)/sqrt(double(j-i)); 
          l_monopole(i) += zdensity(j)/(sqrt(double(j-i))*double(j-i)); // longitudinal monopole contribution
       } 
       //  bunch bunch impedance,  same turn the leading buckets effect
       for (int ibunch= bunch_i+1; ibunch< numbunches; ++ibunch){ 
           double denominator_i =  1.0/
                       (N_factor*sqrt(rbunch_sp*(stored_buckets(ibunch,0)-stored_buckets(bunch_i,0))+num_slices-i-1));        
           
           
           dipole_x(i) +=  stored_bunchnp(ibunch,0)* stored_means(ibunch,0,0)*denominator_i;
           dipole_y(i) +=  stored_bunchnp(ibunch,0)*stored_means(ibunch,0,1)*denominator_i;
           quad(i) += stored_bunchnp(ibunch,0)*denominator_i;

           
           l_monopole(i) +=stored_bunchnp(ibunch,0)/(N_factor*pow(
                      rbunch_sp*(stored_buckets(ibunch,0)-stored_buckets(bunch_i,0))  // distance between buckets 
                   +(stored_means(ibunch,0,2)-stored_means(bunch_i,0,2))*rescaling_f  // relative shift of the  z_means
                   +num_slices-i-1,                                          //distance from slice to the end of the bunch
                      1.5));
      } 
      // bunch-bunch impedance, the previous turn, the following buckets effect
      if (registered_turns>1) {

            for (int ibunch= bunch_i-1; ibunch>= 0; --ibunch){
                double denominator_i =  1.0/
                      (N_factor*sqrt(rbunch_sp*(stored_buckets(ibunch,1)-stored_buckets(bunch_i,0))+line_length*rescaling_f+i-1));
          
                dipole_x(i) +=  stored_bunchnp(ibunch,1)* stored_means(ibunch,1,0)*denominator_i;
                dipole_y(i) +=  stored_bunchnp(ibunch,1)*stored_means(ibunch,1,1)*denominator_i;
                quad(i) += stored_bunchnp(ibunch,1)*denominator_i;

          
                l_monopole(i) +=stored_bunchnp(ibunch,1)/(N_factor*pow(
                     rbunch_sp*(stored_buckets(ibunch,1)-stored_buckets(bunch_i,0))  // distance between buckets 
                  +(stored_means(ibunch,1,2)-stored_means(bunch_i,0,2))*rescaling_f  // relative shift of the  z_means
                 +line_length*rescaling_f
                  +i-1,                                          //distance from slice to the end of the bunch
                  1.5)); 
            } 
      }         
      
    }
//  impedance contribution from previous turns 
    double xdipole=0.;
    double ydipole=0.;
    double quadcontrib=0.;
    double lmonopole=0.;   
 // finishing with previous turn    
// contribution of the leading buckets effect
    if (registered_turns>1) {
        for (int ibunch= bunch_i; ibunch< numbunches; ++ibunch){
            double denominator_i =  1.0/
                    (N_factor*sqrt(rbunch_sp*(stored_buckets(ibunch,1)-stored_buckets(bunch_i,0))+line_length*rescaling_f));
            xdipole += stored_bunchnp(ibunch,1)* stored_means(ibunch,1,0)*denominator_i; 
            ydipole += stored_bunchnp(ibunch,1)* stored_means(ibunch,1,1)*denominator_i;
            quadcontrib += stored_bunchnp(ibunch,1)*denominator_i;
            
            lmonopole +=  stored_bunchnp(ibunch,1)/(N_factor*pow(
                                     rbunch_sp*(stored_buckets(ibunch,1)-stored_buckets(bunch_i,1))  // distance between buckets 
                                             +(stored_means(ibunch,1,2)-stored_means(bunch_i,0,2))*rescaling_f  // relative shift of the  z_means
                                             +line_length*rescaling_f,                      //distance from slice to the end of the bunch
                                             1.5)); 
        }    
    }
// contribution of previous (>1) turns   
    for (int iturn=2; iturn<registered_turns; ++iturn){
        for (int ibunch= 0; ibunch< numbunches; ++ibunch){
            double denominator_i =  1.0/
                    (N_factor*sqrt(rbunch_sp*(stored_buckets(ibunch,iturn)-stored_buckets(bunch_i,0))+iturn*line_length*rescaling_f));
            xdipole += stored_bunchnp(ibunch,iturn)* stored_means(ibunch,iturn,0)*denominator_i; 
            ydipole += stored_bunchnp(ibunch,iturn)* stored_means(ibunch,iturn,1)*denominator_i;
            quadcontrib += stored_bunchnp(ibunch,iturn)*denominator_i;
            
            lmonopole +=  stored_bunchnp(ibunch,iturn)/(N_factor*pow(
                                         rbunch_sp*(stored_buckets(ibunch,iturn)-stored_buckets(bunch_i,iturn))  // distance between buckets 
                                          +(stored_means(ibunch,iturn,2)-stored_means(bunch_i,0,2))*rescaling_f  // relative shift of the  z_means
                                                 +iturn*line_length*rescaling_f,                      //distance from slice to the end of the bunch
                                                 1.5)); 
        } 
    }
        
    dipole_x.add(xdipole);
    dipole_y.add(ydipole);
    quad.add(quadcontrib);
    l_monopole.add(lmonopole);
//         
}

