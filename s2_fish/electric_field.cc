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
    const double mass = mbs.mass * 1.0e9;
    const double eps0 = PH_MKS_eps0; 
    const double qe=PH_MKS_e;
/*    consider the electric fields in the rest frame, E'x, E'y, E'z,
	and Ex,Ey,Ez in the lab frame
      the force on a charge q is
	Fx=q * E'x/gamma =q * Ex/gamma^2
	Fy=q * E'y/gamma=q * Ey/gamma^2
	Fz=q * E'z = q *Ez	
	
	the kick in the 3rd direction is a kick of p_t, not p_z	
	p_t=-U ==> dp_t/dt=-betaz * dpz/dt-betax *dpx/dt-betay*dpy/dt =~ -beta*dpz/dt
*/

    //double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
    // factor *= length* mbs.total_current /(beta * c); // total charge=linear charge density*length
    double factor =PH_CNV_brho_to_p/eps0; // charge of the particle is PH_CNV_brho_to_p =p/Brho in [GeV/c]/[Tesla*m]
    factor *=mbs.bunch_np*qe*mbs.charge;// total charge
    //std::cout<<" previous factor= "<<length* mbs.total_current /(beta * c)<<" now ="<<mbs.bunch_np*qe<<std::endl;
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
    const double qe=PH_MKS_e;

    //double length=2.0*pi*beta/mbs.units(0); // bunch length in lab frame
   
    

   // factor *= length* mbs.total_current /(beta * c); // total charge=linear charge density*length
   // std::cout<<" previous factor= "<<length* mbs.total_current /(beta * c)<<" now ="<<mbs.bunch_np*qe<<std::endl;
    double factor =PH_CNV_brho_to_p/eps0; // charge of the particle is PH_CNV_brho_to_p =p/Brho
    factor *=mbs.bunch_np*qe*mbs.charge;// total charge
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
     const double qe=PH_MKS_e;
     
/*    consider the electric fields in the rest frame, E'x, E'y, E'z,
	and Ex,Ey,Ez in the lab frame
      the force on a charge q is
	Fx=q * E'x/gamma =q * Ex/gamma^2
	Fy=q * E'y/gamma=q * Ey/gamma^2
	Fz=q * E'z = q *Ez	
	
	the kick in the 3rd direction is a kick of p_t, not p_z	
	p_t=-U ==> dp_t/dt=-betaz * dpz/dt-betax *dpx/dt-betay*dpy/dt =~ -beta*dpz/dt
*/

   
    double factor =PH_CNV_brho_to_p/eps0; // charge of the particle is PH_CNV_brho_to_p =p/Brho
    factor *=mbs.bunch_np*qe*mbs.charge;// total charge
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

