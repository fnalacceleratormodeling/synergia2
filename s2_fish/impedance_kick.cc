#include "mpi.h"
#include "impedance_kick.h"
#include <cmath>
#include "math_constants.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "basic_toolkit/PhysicsConstants.h"

void
calculate_rwvars(Macro_bunch_store& mbs,
                 Array_1d<double> &zdensity,
                 Array_1d<double> &xmom, Array_1d<double> &ymom,
                 double z_left, double z_length, Array_1d<int> &bin_partition)
{


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
  //   this routine is expecting particles in the fixed-t frame with
  //   the units factor already applied
//  const double really_big = 1.0e30;
//  double zmin = really_big;
//  double zmax = -really_big;
//    for (int n = 0; n < mbs.local_num; ++n) {
//      double z = mbs.local_particles(4,n);
//      if (z < zmin) {
//          zmin = z;
//      } else if (z > zmax) {
//          zmax = z;
//      }
//  }
//    double z_left = zmin;
//    double z_length = (zmax - zmin);
//    std::cout << "jfa: z_left = " << z_left << ", z_length = " << z_length << std::endl;
    int z_num = zdensity.get_length();
//    double h = z_length/(z_num-1.0); // AM , before was h = z_length/z_num
    double h = z_length/z_num;  // AM, what if z_periodic??? 
    Array_1d<double> local_zdensity(z_num);
    Array_1d<double> local_xmom(z_num);
    Array_1d<double> local_ymom(z_num);

    local_zdensity.set_all(0.0);
    local_xmom.set_all(0.0);
    local_ymom.set_all(0.0);
    for (int n = 0; n < mbs.local_num; ++n) {
        int bin = static_cast<int>((mbs.local_particles(4,n)-z_left)/h);
        if ((bin < z_num) && (bin >= 0)) {
            local_zdensity(bin) += 1;
            local_xmom(bin) += mbs.local_particles(0,n);
            local_ymom(bin) += mbs.local_particles(2,n);
            bin_partition(n)=bin; //bin_partition(n) is the bin where you find the particle n
        }
        else if ((bin==z_num) && fabs(mbs.local_particles(4,n)-z_length-z_left)<z_length*1.e-14) { 
            local_zdensity(bin-1) += 1; // put the edge particle in the last bin=z_num-1
            bin_partition(n)=z_num-1;
        } 
        else
        {   std::cout << "  z_left  "<<z_left<<"  rank= "<<rank<<std::endl;
        std::cout<<"mbs.local_particles(4,n)="  <<mbs.local_particles(4,n)<<"  rank= "<<rank<<std::endl; 
        std::cout<< "NNNNNNNNNN n="<<n<<std::endl; 
        std::cout << "  z_length  "<<z_length<<"  rank= "<<rank<<std::endl;
        std::cout << "(mbs.local_particles(4,n)-z_left)= "<<(mbs.local_particles(4,n)-z_left)<<"  rank= "<<rank<<std::endl;
        std::cout << "bin: " << bin<<"  z_num="<<z_num<< "  h=" << h <<"  rank= "<<rank<<std::endl;                
        std::cout << "mbs.local_particles(4,n)-z_length-z_left= "<<fabs(mbs.local_particles(4,n)-z_length-z_left)<<"  rank= "<<rank<<std::endl;
        throw
                std::runtime_error("particles out of range in calculate_rwvars ");
        }

    
    }
    
    // jfa: the other deposit functions do not do the communication. This one does.
    // Obviously, something should change.
    MPI_Allreduce(reinterpret_cast<void*>(local_zdensity.get_data_ptr()),
                  reinterpret_cast<void*>(zdensity.get_data_ptr()),
                                          z_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(reinterpret_cast<void*>(local_xmom.get_data_ptr()),
                  reinterpret_cast<void*>(xmom.get_data_ptr()),
                                          z_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(reinterpret_cast<void*>(local_ymom.get_data_ptr()),
                  reinterpret_cast<void*>(ymom.get_data_ptr()),
                                          z_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //    dbg is set here XXXXXXXXXXXXXX
    int dbg = 0;

    for (int k = 0; k < z_num; ++k) {
        if (zdensity(k) != 0.0) {
            if (dbg) std::cout << "before bin: " << k << " zdensity(k): " << zdensity(k) << " xmom(k): " << xmom(k) << " ymom(k): " << ymom(k) << " units: " << mbs.units(0) << std::endl;
            xmom(k) /= zdensity(k);
            ymom(k) /= zdensity(k);
            if (dbg) std::cout << "after bin: " << k << " zdensity(k): " << zdensity(k) << " xmom(k): " << xmom(k) << " ymom(k): " << ymom(k) << std::endl;
        } else {
            xmom(k) = 0.0;
            ymom(k) = 0.0;
        }
    }
    
}   



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
        int registered_turns=stored_means.get_shape()[1];
         try{
             if (fabs((stored_bunchnp(bunch_i,0)-stored_bunchnp(bunch_i,registered_turns-1)))>stored_bunchnp(bunch_i,0)*1.e-16){
                 throw std::runtime_error("the bunchnp changes with turns, quad impedance assumes it is constant!");
             }
         } catch ( std::runtime_error e){
             std::clog<<e.what()<<std::endl;
        }

        double quad_wake_sum_scaled=stored_bunchnp(bunch_i,registered_turns-1)*quad_wake_sum/(N_factor*sqrt(rescaling_f)); //rescaled to cancel
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

