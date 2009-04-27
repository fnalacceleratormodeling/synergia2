#include "solver_fftw.h"
#include <iomanip>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

#undef DL_IMPORT
#include "mytimer.h"
#include "fftw_helper.h"

#include "communicate.h"

Complex_scalar_field
get_rho_hat2(Real_scalar_field &rho, Fftw_helper &fftwh)
{
    // steps 1 and 2
    Int3 num_points2 = rho.get_points().get_shape();
   // num_points2.scale(2); am: do we need this?
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field rho2(fftwh.padded_shape_real().vector(),
                           physical_size2.vector(),
                           rho.get_physical_offset(),
                           fftwh.guard_lower(), fftwh.guard_upper());
    //  rho2.get_points().set_storage_size(fftwh.local_size());
    timer("misc");
    rho2.get_points().zero_all();
    timer("calc zero all rho2");
    Int3 num_points = rho.get_points().get_shape();
    Int3 index;
    timer("misc");
    int index0_max = std::min(fftwh.upper(), num_points[0]);
    
    for (index[0] = fftwh.lower(); index[0] < index0_max; ++index[0]) {
        for (index[1] = 0; index[1] < num_points[1]; ++index[1]) {
            for (index[2] = 0; index[2] < num_points[2]; ++index[2]) {
                rho2.get_points().set(index, rho.get_points().get(index));
            }
        }
    }
    timer("calc rho2");

    Complex_scalar_field rho_hat2(fftwh.padded_shape_complex().vector(),
                                  physical_size2.vector(),
                                  rho.get_physical_offset(),
                                  fftwh.guard_lower(), fftwh.guard_upper());
    //  rho_hat2.get_points().set_storage_size(fftwh.local_size());
    fftwh.transform(rho2, rho_hat2);
    timer("rho fft");
    return rho_hat2;
}

Real_scalar_field
get_G2(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh)
{
    const double pi = 4.0 * atan(1.0);
    Int3 num_points = rho.get_points().get_shape();
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
    Double3 physical_size = rho.get_physical_size();
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field G2(fftwh.padded_shape_real().vector(),
                         physical_size2.vector(),
                         rho.get_physical_offset(),
                         fftwh.guard_lower(), fftwh.guard_upper());
    //  G2.get_points().set_storage_size(fftwh.local_size());
    Double3 h(rho.get_cell_size());
    Int3 index;
    double num_points_z;	
    (!z_periodic) ? num_points_z=num_points[2] : num_points_z=num_points[2]-1;

    // What is G(0,0,0), anyway? Rob and Ji seem to think it is G(0,0,1).
    // Hockney seems to think it is 1.
    // I don't think it is either, but I have not yet worked out what I
    // consider to be the right answer (the one that preserves the integral
    // of G*rho).

    // Rob and Ji version
    //   G2.get_points().set(Int3(0,0,0),G2.get_points().get(Int3(0,0,1)));

    // This would be the correct value if we were using cells that were spheres
    // instead of rectangular solids:
    // average value of inner sphere:
    //   G2.get_points().set(Int3(0,0,0),(1.0/4.0*pi)*
    // 		      (3.0/(2.0*((1+sqrt(3))/2.0)*
    // 			    sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]))));

    // average value of outer sphere. This works unreasonably well.
    double G000 = (1.0 / 4.0 * pi) * (3.0 / (2.0 * (sqrt(3)) *
                                      sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2])));
    // Calculate what I think the answer should be by doing a very
    // simple numerical integral. The resulting answer doesn't work
    // nearly as well as the outer sphere approximation above, so
    // something must be wrong.
    //   double sum = 0.0;
    //   const int nsteps = 32;
    //   double dx = h[0]/nsteps;
    //   double dy = h[1]/nsteps;
    //   double dz = h[2]/nsteps;
    //   x = 0.5*dx;
    //   for(int i=0; i<nsteps; ++i) {
    //     y = 0.5*dy;
    //     for(int j=0; j<nsteps; ++j) {
    //       z = 0.5*dz;
    //       for(int k=0; k<nsteps; ++k) {
    // 	sum += 1.0/sqrt(x*x+y*y+z*z);
    // 	z += dz;
    //       }
    //       y += dy;
    //     }
    //     x += dx;
    //   }
    // //   double R = h[0];
    // //   double vol = (1/8.0)*4.0*pi/3.0*R*R*R;
    //   double R = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    //   double vol = h[0]*h[1]*h[2];
    //   double mean_inv_r = sum*dx*dy*dz/vol;
    //   G2.get_points().set(Int3(0,0,0),(1.0/4.0*pi)*mean_inv_r);
    timer("misc");
    double x, y, z, G;
    const int num_images = 8;
    int miy, miz; // mirror index x, etc.
   // double z_bin_offset;
   // z_bin_offset = 0.0;
//     if (z_periodic) {
//         z_bin_offset = 0.0;//-0.5;
//     } else {
//         z_bin_offset = 0.0;
//     }


    for (index[0] = fftwh.lower(); index[0] < fftwh.upper(); ++index[0]) {
        if (index[0] > num_points2[0] / 2) {
            x = (num_points2[0] - index[0]) * h[0];
        } else {
            x = index[0] * h[0];
        }
        for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
            y = index[1] * h[1];
            miy = num_points2[1] - index[1];
            if (miy == num_points[1]) {
                miy = num_points2[1];
            }
            for (index[2] = 0; index[2] <= num_points_z; ++index[2]) {

                 z = index[2]* h[2];            
	         if (!z_periodic) {miz=num_points2[2]-index[2]; 
                                if (miz == num_points[2])  miz = num_points2[2];
                  }              
                if (!((x == 0.0) && (y == 0.0) && (z == 0.0))) {
                    G = 1.0 / (4.0 * pi * sqrt(x * x + y * y + z * z));
                } else {
                    G = G000;
                }
                if (z_periodic) {
                    for (int image = -num_images; image < num_images; ++image) {
                        if (image != 0) {
                            double z_image = z + image * physical_size[2];
                            if (!((x == 0.0) && (y == 0.0) && (fabs(z_image) < 1.0e-14))) {
                                G += 1.0 / (4.0 * pi * sqrt(x * x + y * y + z_image * z_image));
                            } else {
                                G += G000;
                            }
                        }
                    }
                }
                G2.get_points().set(index, G);
                // three mirror images
                if (miy < num_points2[1]) {
                    G2.get_points().set(Int3(index[0], miy, index[2]), G);
                    if ((miz < num_points2[2]) && (!z_periodic)){
                        G2.get_points().set(Int3(index[0], miy, miz), G);
                    }
                }
                if ((miz < num_points2[2]) && (!z_periodic)){
                    G2.get_points().set(Int3(index[0], index[1], miz), G);
                }
            }
        }
    }
    timer("calc G");
    return G2;
}

Real_scalar_field
get_G2_z_steps(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh)
{

	

    const double pi = 4.0 * atan(1.0);
    int distribute_fftwh;
    Int3 num_points = rho.get_points().get_shape();
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
//      if ((fftwh.lower()!= 0) ||(fftwh.upper() != num_points2[0])){
//         throw
//                std::runtime_error("get_G2_z_steps requires the guards fftwh.lower()= 0 and fftwh.upper()=num_points2[0]), and this is not the case");
//     }

    int x_lrange=fftwh.lower();
    int x_urange=fftwh.upper();
    if ((x_lrange != 0) ||(x_urange != num_points2[0])) { 
				distribute_fftwh=1;     }
    else {                      distribute_fftwh=0;
                                x_urange=x_urange/2+1;     }

    double num_points_z;	
    (!z_periodic) ? num_points_z=num_points[2] : num_points_z=num_points[2]-1;

    Double3 physical_size = rho.get_physical_size();
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field G2(fftwh.padded_shape_real().vector(),
                         physical_size2.vector(),
                         rho.get_physical_offset(),
                         fftwh.guard_lower(), fftwh.guard_upper());
    //  G2.get_points().set_storage_size(fftwh.local_size());
    Double3 h(rho.get_cell_size());
    Int3 index;
   
    //r1=sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
   // double G000 = (2.0* pi)*(h[0] * h[0] + h[1] * h[1] + h[2] * h[2])/(4*pi*r1*r1*r1);


    timer("misc");
    double x, y, z, G, G000;
    const int num_images = 8;
    int mix, miy, miz; // mirror index x, etc.
    double z_bin_offset, hz, rr, r1,r2;
    hz=1.*h[2];
    rr= h[0] * h[0] + h[1] * h[1];
    r1=sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
    G000=(2.0/rr)*(h[2]*r1+rr*log((h[2]+r1)/sqrt(rr))-h[2] * h[2]);// average value of outer cylinder.

        G=0.;
	

	
    //    for (index[0] = 0; index[0] <= num_points2[0]/2; ++index[0]) {
         for (index[0] = x_lrange; index[0] <x_urange; ++index[0]) {
            if (index[0] > num_points2[0] / 2) {
                x = (num_points2[0] - index[0]) * h[0];
            } else {
                 x = index[0] * h[0];
            }             
            if (distribute_fftwh==0) {mix=num_points2[0]-index[0];
                                     if (mix == num_points[0])  mix = num_points2[0];}
            
            for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
	        y = index[1] * h[1];
	        miy=num_points2[1]-index[1];
                rr=x * x + y * y;
	        if (miy == num_points[1])  miy = num_points2[1];

                for (index[2] = 0; index[2] <= num_points_z; ++index[2]) {
	           z = index[2]* h[2];            
	           if (!z_periodic) {miz=num_points2[2]-index[2]; 
                                if (miz == num_points[2])  miz = num_points2[2];
                   }


             if (fabs(fabs(z)-hz)<1.0e-10){   // z=+-hz  case
                        if (fabs(x)+fabs(y)<2.0e-10) { // r_perp=0
				r1=(sqrt(rr+(fabs(z)+hz)*(fabs(z)+hz))+fabs(z)+hz)/(sqrt(rr+h[2]*h[2])+h[2]);
				G=G000/2.+log(r1);   }
                        else {r1=(sqrt(rr+(fabs(z)+hz)*(fabs(z)+hz))+fabs(z)+hz)/sqrt(rr);// r_perp != 0
                             G=log(r1);}
                        }
             else if ((z<hz-1.0e-10)&&(z>-hz+1.0e-10)) { //-hz<z<hz  case		
	        if (fabs(x)+fabs(y)<2.0e-10) {G=G000;}   // r_perp=0  case
		else {r1=(sqrt(rr+h[2]*h[2])+h[2])/sqrt(rr); // r_perp != 0
		      G=2.*log(r1);}
		r1=(sqrt(rr+(z-hz)*(z-hz))-z+hz)/(sqrt(rr+h[2]*h[2])+h[2]);
		r2=(sqrt(rr+(z+hz)*(z+hz))+z+hz)/(sqrt(rr+h[2]*h[2])+h[2]);
		G += log(r1)+log(r2);	
                }
            else { //  z<-hz  or hz<z case
                     r1=(sqrt((fabs(z)+hz)*(fabs(z)+hz)+rr)+fabs(z)+hz)/
                      (sqrt((fabs(z)-hz)*(fabs(z)-hz)+rr)+fabs(z)-hz);
                      G = log(fabs(r1));
                                         }


             if (z_periodic) {
               for (int image = -num_images; image < num_images; ++image) {
                   if (image != 0) {
                   double z_image = z + image * physical_size[2];

                    if (fabs(fabs(z_image)-hz)<1.0e-10){
                        if (fabs(x)+fabs(y)<2.0e-10) {
				r1=(sqrt(rr+(fabs(z_image)+hz)*(fabs(z_image)+hz))+fabs(z_image)+hz)/(sqrt(rr+h[2]*h[2])+h[2]);
				G +=G000/2.+log(r1);   }
                        else {r1=(sqrt(rr+(fabs(z_image)+hz)*(fabs(z_image)+hz))+fabs(z_image)+hz)/sqrt(rr);
                             G +=log(r1);}
                        }
                    else if ((z_image<hz-1.0e-10)&&(z_image>-hz+1.0e-10)) {		
	               if (fabs(x)+fabs(y)<2.0e-10) {G=+G000;}
		       else {r1=(sqrt(rr+h[2]*h[2])+h[2])/sqrt(rr);
		             G +=2.*log(r1);}
		        r1=(sqrt(rr+(z_image-hz)*(z_image-hz))-z_image+hz)/(sqrt(rr+h[2]*h[2])+h[2]);
		        r2=(sqrt(rr+(z_image+hz)*(z_image+hz))+z_image+hz)/(sqrt(rr+h[2]*h[2])+h[2]);
		        G += log(r1)+log(r2);	
                             }
                   else {
                          r1=(sqrt((fabs(z_image)+hz)*(fabs(z_image)+hz)+rr)+fabs(z_image)+hz)/
                          (sqrt((fabs(z_image)-hz)*(fabs(z_image)-hz)+rr)+fabs(z_image)-hz);
                           G += log(fabs(r1));
                                         }


                       }
                    }
                }
		
                G2.get_points().set(index, G);



                // seven mirror images
                if ((distribute_fftwh==0)  && (mix < num_points2[0])) {
                    G2.get_points().set(Int3(mix, index[1], index[2]), G);
                    if (miy < num_points2[1]) {
                        G2.get_points().set(Int3(mix, miy, index[2]), G);
                        if ((!z_periodic)  &&(miz < num_points2[2])) {
                            G2.get_points().set(Int3(mix, miy, miz), G);
                        }
                    }
                    if ((!z_periodic) && (miz < num_points2[2])) {
                        G2.get_points().set(Int3(mix, index[1], miz), G);
                    }
                }
                if (miy < num_points2[1]) {
                    G2.get_points().set(Int3(index[0], miy, index[2]), G);
                    if ((!z_periodic) && (miz < num_points2[2])) {
                        G2.get_points().set(Int3(index[0], miy, miz), G);
                    }
                }
                if ((!z_periodic) && (miz < num_points2[2])) {
                    G2.get_points().set(Int3(index[0], index[1], miz), G);
                }

            }
        }
    }
    double scale=1.0/(4.0 * pi*2.0*hz);	
    G2.get_points().scale(scale);
    timer("calc G");
    return G2;
}

Real_scalar_field
get_G2_z_linear(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh)
{

	

    const double pi = 4.0 * atan(1.0);
    int distribute_fftwh;
    Int3 num_points = rho.get_points().get_shape();
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
    

    int x_lrange=fftwh.lower();
    int x_urange=fftwh.upper();
    if ((x_lrange != 0) ||(x_urange != num_points2[0])) { 
				distribute_fftwh=1;     }
    else {                      distribute_fftwh=0;
                                x_urange=x_urange/2+1;     }
    double num_points_z;	
    (!z_periodic) ? num_points_z=num_points[2] : num_points_z=num_points[2]-1;

    Double3 physical_size = rho.get_physical_size();
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field G2(fftwh.padded_shape_real().vector(),
                         physical_size2.vector(),
                         rho.get_physical_offset(),
                         fftwh.guard_lower(), fftwh.guard_upper());
    //  G2.get_points().set_storage_size(fftwh.local_size());
    Double3 h(rho.get_cell_size());
    Int3 index;
   


    timer("misc");
    double x, y, z, G, G000;
    const int num_images = 8;
    int mix, miy, miz; // mirror index x, etc.
    double z_bin_offset, rr, r1,r2, T1,T2;
    const double hz=h[2]; // do not modify hz!
	


     rr= h[0] * h[0] + h[1] * h[1];
     r1=sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
     G000=(2.0/rr)*(h[2]*r1+rr*log((h[2]+r1)/sqrt(rr))-h[2] * h[2]);// average value of outer cylinder.

     
     G=0.;
      // for (index[0] = 0; index[0] <= num_points2[0]/2; ++index[0]) {
        for (index[0] = x_lrange; index[0] <x_urange; ++index[0]) {
           if (index[0] > num_points2[0] / 2) {
                x = (num_points2[0] - index[0]) * h[0];
           } else {
                 x = index[0] * h[0];
           }        
           if (distribute_fftwh==0) {mix=num_points2[0]-index[0];
                                     if (mix == num_points[0])  mix = num_points2[0];
           }
            
           for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
	     y = index[1] * h[1];
	     miy=num_points2[1]-index[1];
	     if (miy == num_points[1])  miy = num_points2[1];
             rr=x * x + y * y;
             for (index[2] = 0; index[2] <= num_points_z; ++index[2]) {
	         z = index[2]* h[2];            
	         if (!z_periodic) {miz=num_points2[2]-index[2]; 
                                if (miz == num_points[2])  miz = num_points2[2];
                 }

	        
		 G=2.0*sqrt(rr+z*z)-sqrt(rr+(z-hz)*(z-hz))-sqrt(rr+(z+hz)*(z+hz));
		 if (z<-hz) {
                           r1=(sqrt((z-hz)*(z-hz)+rr)-z+hz)/(sqrt(z*z+rr)-z);
			   T1=(hz-z)*log(r1);
			   r2=(sqrt(z*z+rr)-z)/(sqrt((z+hz)*(z+hz)+rr)-z-hz);
			   T2=(hz+z)*log(r2);
			   G +=	T1+T2; }

		else if (fabs(z+hz)<1.0e-10){
			   r1=(sqrt((z-hz)*(z-hz)+rr)-z+hz)/(sqrt(z*z+rr)-z);
			   T1=(hz-z)*log(r1);
			   G +=	T1;}
		
 		else if (fabs(z)<1.e-10){
			   if (fabs(x)+fabs(y)<2.0e-10)	{G +=hz*G000;} // T1+T2 in fact		
			   else {                        r1=(sqrt(hz*hz+rr)+hz)/sqrt(rr);
				                         G += 2.*hz*log(r1);}	
					}
		else if (fabs(z-hz)<1.0e-10){ 
			   r1=(sqrt((z+hz)*(z+hz)+rr)+z+hz)/(sqrt(z*z+rr)+z);
			   T1=(hz+z)*log(r1);
			   G +=	T1;}
		else if (z>hz){
			   r1=(sqrt(z*z+rr)+z)/(sqrt((z-hz)*(z-hz)+rr)+z-hz);			   
			   T1=(hz-z)*log(r1);
                           r2=(sqrt((z+hz)*(z+hz)+rr)+z+hz)/(sqrt(z*z+rr)+z);			
			   T2=(hz+z)*log(r2);
			   G +=	T1+T2;}
		else{ throw
                    std::runtime_error("get_G2_z_linear error1, check if hz=h[2]");}

             
                if (z_periodic) {
                   for (int image = -num_images; image < num_images; ++image) {
                        if (image != 0) {
                           double z_image = z + image * physical_size[2];


			if (z_image<-hz) {
                           r1=(sqrt((z_image-hz)*(z_image-hz)+rr)-z_image+hz)/(sqrt(z_image*z_image+rr)-z_image);
			   T1=(hz-z_image)*log(r1);
			   r2=(sqrt(z_image*z_image+rr)-z_image)/(sqrt((z_image+hz)*(z_image+hz)+rr)-z_image-hz);
			   T2=(hz+z_image)*log(r2);
			   G +=	T1+T2; }

			else if (fabs(z_image+hz)<1.0e-10){
			   r1=(sqrt((z_image-hz)*(z_image-hz)+rr)-z_image+hz)/(sqrt(z_image*z_image+rr)-z_image);
			   T1=(hz-z_image)*log(r1);
			   G +=	T1;}
		
 			else if (fabs(z_image)<1.e-10){
			   if (fabs(x)+fabs(y)<2.0e-10)	{G +=hz*G000;} // T1+T2 in fact		
			   else {                        r1=(sqrt(hz*hz+rr)+hz)/sqrt(rr);
				                         G += 2.*hz*log(r1);}	
					}
			else if (fabs(z_image-hz)<1.0e-10){ 
			   r1=(sqrt((z_image+hz)*(z_image+hz)+rr)+z_image+hz)/(sqrt(z_image*z_image+rr)+z_image);
			   T1=(hz+z_image)*log(r1);
			   G +=	T1;}
			else if (z_image>hz){
			   r1=(sqrt(z_image*z_image+rr)+z_image)/
					(sqrt((z_image-hz)*(z_image-hz)+rr)+z_image-hz);			   
			   T1=(hz-z_image)*log(r1);
                           r2=(sqrt((z_image+hz)*(z_image+hz)+rr)+z_image+hz)/
							(sqrt(z_image*z_image+rr)+z_image);			
			   T2=(hz+z_image)*log(r2);
			   G +=	T1+T2;}
			else{ throw
                               std::runtime_error("get_G2_z_linear error2, check if hz=h[2]");}
   			 

                       }
                    }
                }
                G2.get_points().set(index, G);


                // seven mirror images
                if ((distribute_fftwh==0)  && (mix < num_points2[0])) {
                    G2.get_points().set(Int3(mix, index[1], index[2]), G);
                    if (miy < num_points2[1]) {
                        G2.get_points().set(Int3(mix, miy, index[2]), G);
                        if ((!z_periodic)  &&(miz < num_points2[2])) {
                            G2.get_points().set(Int3(mix, miy, miz), G);
                        }
                    }
                    if ((!z_periodic) && (miz < num_points2[2])) {
                        G2.get_points().set(Int3(mix, index[1], miz), G);
                    }
                }
                if (miy < num_points2[1]) {
                    G2.get_points().set(Int3(index[0], miy, index[2]), G);
                    if ((!z_periodic) && (miz < num_points2[2])) {
                        G2.get_points().set(Int3(index[0], miy, miz), G);
                    }
                }
                if ((!z_periodic) && (miz < num_points2[2])) {
                    G2.get_points().set(Int3(index[0], index[1], miz), G);
                }


            }
        }
    }
    double scale=1.0/(4.0 * pi*hz*hz);	
    G2.get_points().scale(scale);
    timer("calc G");
    return G2;
}

Real_scalar_field
get_G2_spherical(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh)
{

	

    const double pi = 4.0 * atan(1.0);
    int distribute_fftwh;
    Int3 num_points = rho.get_points().get_shape();
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
   
     
    int x_lrange=fftwh.lower();
    int x_urange=fftwh.upper();
    if ((x_lrange != 0) ||(x_urange != num_points2[0])) { 
				distribute_fftwh=1;     }
    else {                      distribute_fftwh=0;
                                x_urange=x_urange/2+1;     }


    double num_points_z;	
    (!z_periodic) ? num_points_z=num_points[2] : num_points_z=num_points[2]-1;

	
    Double3 physical_size = rho.get_physical_size();
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
  // if (z_periodic) { physical_size2[2] /= 2.0;}
    Real_scalar_field G2(fftwh.padded_shape_real().vector(),
                         physical_size2.vector(),
                         rho.get_physical_offset(),
                         fftwh.guard_lower(), fftwh.guard_upper());
    //  G2.get_points().set_storage_size(fftwh.local_size());
    Double3 h(rho.get_cell_size());
    Int3 index;
   


    timer("misc");
    double x, y, z, rr, G, G0,G1;
    const int num_images = 8;
    int mix, miy, miz; // mirror index x, etc.
    double hr=0.5*sqrt(h[0] * h[0] + h[1]*h[1]);//Note that h[2] usually much larger than h[0] and h[1],
						// so probably this is not going to work too well!

     G0=hr*hr/2.0; 
     G1=hr*hr*hr/3.0;	   
     G=0.;
 
        //for (index[0] = 0; index[0] <= num_points[0]; ++index[0]) {
         for (index[0] = x_lrange; index[0] <x_urange; ++index[0]) {          
             if (index[0] > num_points2[0] / 2) {
                x = (num_points2[0] - index[0]) * h[0];
             } else {
                 x = index[0] * h[0];
             }
             if (distribute_fftwh==0) {mix=num_points2[0]-index[0];
                                     if (mix == num_points[0])  mix = num_points2[0];
             }
	  
            
           for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
	     y = index[1] * h[1];
	     miy=num_points2[1]-index[1];
	     if (miy == num_points[1])  miy = num_points2[1];
             for (index[2] = 0; index[2] <= num_points_z; ++index[2]) {
	       z = index[2]* h[2];            
	       if (!z_periodic) {miz=num_points2[2]-index[2]; 
                                if (miz == num_points[2])  miz = num_points2[2];}

	       
              rr=sqrt(x*x+y*y+z*z);
	      if(rr<1.e-10) {G=G0;}
	      else if (rr<hr){
	      G=G0-rr*rr/6.0;}
              else {G=G1/rr;}
		


             
                if (z_periodic) {
                   for (int image = -num_images; image < num_images; ++image) {
                        if (image != 0) {
                           double z_image = z + image * physical_size[2];

 	      rr=sqrt(x*x+y*y+z_image*z_image);
	      if(rr<1.e-10) {G += G0;}
	      else if (rr<hr){		
              G += G0-rr*rr/6.0;}
              else {G += G1/rr;}


                       }
                    }
                }
                G2.get_points().set(index, G);
		

                // seven mirror images
                if ((distribute_fftwh==0)  && (mix < num_points2[0])) {
                    G2.get_points().set(Int3(mix, index[1], index[2]), G);
                    if (miy < num_points2[1]) {
                        G2.get_points().set(Int3(mix, miy, index[2]), G);
                        if ((!z_periodic)  &&(miz < num_points2[2])) {
                            G2.get_points().set(Int3(mix, miy, miz), G);
                        }
                    }
                    if ((!z_periodic) && (miz < num_points2[2])) {
                        G2.get_points().set(Int3(mix, index[1], miz), G);
                    }
                }
                if (miy < num_points2[1]) {
                    G2.get_points().set(Int3(index[0], miy, index[2]), G);
                    if ((!z_periodic) && (miz < num_points2[2])) {
                        G2.get_points().set(Int3(index[0], miy, miz), G);
                    }
                }
                if ((!z_periodic) && (miz < num_points2[2])) {
                    G2.get_points().set(Int3(index[0], index[1], miz), G);
                }


            }
        }
    }
    double scale=3.0/(4.0*pi*hr*hr*hr);	
    G2.get_points().scale(scale);      
    timer("calc G");
    return G2;
}


Real_scalar_field
get_G2_old(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh)
{
    const double pi = 4.0 * atan(1.0);
    Int3 num_points = rho.get_points().get_shape();
    Int3 num_points2 = rho.get_points().get_shape();
    num_points2.scale(2);
    Double3 physical_size = rho.get_physical_size();
    Double3 physical_size2 = rho.get_physical_size();
    physical_size2.scale(2.0);
    Real_scalar_field G2(fftwh.padded_shape_real().vector(),
                         physical_size2.vector(),
                         rho.get_physical_offset(),
                         fftwh.guard_lower(), fftwh.guard_upper());
    //  G2.get_points().set_storage_size(fftwh.local_size());
    Double3 h(rho.get_cell_size());
    Int3 index;
    // What is G(0,0,0), anyway? Rob and Ji seem to think it is G(0,0,1).
    // Hockney seems to think it is 1.
    // I don't think it is either, but I have not yet worked out what I
    // consider to be the right answer (the one that preserves the integral
    // of G*rho).

    // Rob and Ji version
    //   G2.get_points().set(Int3(0,0,0),G2.get_points().get(Int3(0,0,1)));

    // This would be the correct value if we were using cells that were spheres
    // instead of rectangular solids:
    // average value of inner sphere:
    //   G2.get_points().set(Int3(0,0,0),(1.0/4.0*pi)*
    // 		      (3.0/(2.0*((1+sqrt(3))/2.0)*
    // 			    sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]))));

    // average value of outer sphere. This works unreasonably well.
    double G000 = (1.0 / 4.0 * pi) * (3.0 / (2.0 * (sqrt(3)) *
                                     sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2])));
    // Calculate what I think the answer should be by doing a very
    // simple numerical integral. The resulting answer doesn't work
    // nearly as well as the outer sphere approximation above, so
    // something must be wrong.
    //   double sum = 0.0;
    //   const int nsteps = 32;
    //   double dx = h[0]/nsteps;
    //   double dy = h[1]/nsteps;
    //   double dz = h[2]/nsteps;
    //   x = 0.5*dx;
    //   for(int i=0; i<nsteps; ++i) {
    //     y = 0.5*dy;
    //     for(int j=0; j<nsteps; ++j) {
    //       z = 0.5*dz;
    //       for(int k=0; k<nsteps; ++k) {
    // 	sum += 1.0/sqrt(x*x+y*y+z*z);
    // 	z += dz;
    //       }
    //       y += dy;
    //     }
    //     x += dx;
    //   }
    // //   double R = h[0];
    // //   double vol = (1/8.0)*4.0*pi/3.0*R*R*R;
    //   double R = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    //   double vol = h[0]*h[1]*h[2];
    //   double mean_inv_r = sum*dx*dy*dz/vol;
    //   G2.get_points().set(Int3(0,0,0),(1.0/4.0*pi)*mean_inv_r);
    timer("misc");
    double x, y, z, G;
    const int num_images = 8;
    int mix, miy, miz; // mirror index x, etc.
    for (index[0] = 0; index[0] <= num_points[0]; ++index[0]) {
        x = index[0] * h[0];
        mix = num_points2[0] - index[0];
        if (mix == num_points[0]) {
            mix = num_points2[0];
        }
        for (index[1] = 0; index[1] <= num_points[1]; ++index[1]) {
            y = index[1] * h[1];
            miy = num_points2[1] - index[1];
            if (miy == num_points[1]) {
                miy = num_points2[1];
            }
            for (index[2] = 0; index[2] <= num_points[2]; ++index[2]) {
                z = index[2] * h[2];
                miz = num_points2[2] - index[2];
                if (miz == num_points[2]) {
                    miz = num_points2[2];
                }
                if (!((x == 0.0) && (y == 0.0) && (z == 0.0))) {
                    G = 1.0 / (4.0 * pi * sqrt(x * x + y * y + z * z));
                } else {
                    G = G000;
                }
                if (z_periodic) {
                    for (int image = -num_images; image < num_images; ++image) {
                        if (image != 0) {
                            double z_image = z + image * physical_size[2];
                            if (!((x == 0.0) && (y == 0.0) && (fabs(z_image) < 1.0e-14))) {
                                G += 1.0 / (4.0 * pi * sqrt(x * x + y * y + z_image * z_image));
                            } else {
                                G += G000;
                            }
                        }
                    }
                }
                G2.get_points().set(index, G);
                // seven mirror images
                if (mix < num_points2[0]) {
                    G2.get_points().set(Int3(mix, index[1], index[2]), G);
                    if (miy < num_points2[1]) {
                        G2.get_points().set(Int3(mix, miy, index[2]), G);
                        if (miz < num_points2[2]) {
                            G2.get_points().set(Int3(mix, miy, miz), G);
                        }
                    }
                    if (miz < num_points2[2]) {
                        G2.get_points().set(Int3(mix, index[1], miz), G);
                    }
                }
                if (miy < num_points2[1]) {
                    G2.get_points().set(Int3(index[0], miy, index[2]), G);
                    if (miz < num_points2[2]) {
                        G2.get_points().set(Int3(index[0], miy, miz), G);
                    }
                }
                if (miz < num_points2[2]) {
                    G2.get_points().set(Int3(index[0], index[1], miz), G);
                }
            }
        }
    }
    timer("calc G");
    return G2;
}

Complex_scalar_field
get_G_hat2(Real_scalar_field &rho, bool z_periodic, Fftw_helper &fftwh)
{
    //step 3

     Real_scalar_field G2 = get_G2(rho, z_periodic, fftwh); // AM: after correcting a bug in the z-dimension  
							    //  of the FT, this seems to work the best!
     // Real_scalar_field G2 = get_G2_spherical(rho, z_periodic, fftwh);
    

   // Real_scalar_field G2 = get_G2_z_linear(rho, z_periodic, fftwh);
    // Real_scalar_field G2 = get_G2_z_steps(rho, z_periodic, fftwh);
   
  //  Real_scalar_field G2 = get_G2_old(rho, z_periodic, fftwh);
      Complex_scalar_field G_hat2(fftwh.padded_shape_complex().vector(),
                                G2.get_physical_size(),
                                rho.get_physical_offset(),
                                fftwh.guard_lower(), fftwh.guard_upper());
    //  G_hat2.get_points().set_storage_size(fftwh.local_size());
    fftwh.transform(G2, G_hat2);
    timer("G fft");
    return G_hat2;
}

Complex_scalar_field
get_phi_hat2(Real_scalar_field &rho, Complex_scalar_field &rho_hat2,
             Complex_scalar_field &G_hat2, Fftw_helper &fftwh)
{
    // step 4
    Complex_scalar_field phi_hat2(fftwh.padded_shape_complex().vector(),
                                  G_hat2.get_physical_size(),
                                  G_hat2.get_physical_offset(),
                                  fftwh.guard_lower(), fftwh.guard_upper());
    //  phi_hat2.get_points().set_storage_size(fftwh.local_size());
    Int3 shape(G_hat2.get_points().get_shape());
    Double3 h(rho.get_cell_size());
    timer("misc");

    for (int i = 0; i < G_hat2.get_points().get_length(); ++i) {
        phi_hat2.get_points().get_base_address()[i] =
            rho_hat2.get_points().get_base_address()[i] *
            G_hat2.get_points().get_base_address()[i] *
            h[0] * h[1] * h[2];
    }

    
    timer("calc phi_hat");
    return phi_hat2;
}

Real_scalar_field
get_phi2(Real_scalar_field &rho, Complex_scalar_field &phi_hat2,
         Fftw_helper &fftwh, bool z_periodic)
{
    // step 5
    Int3 num_points2(rho.get_points().get_shape());
    num_points2[0] *= 2;
    num_points2[1] *= 2;
    if (! z_periodic) {
        num_points2[2] *= 2;
    } else {
	num_points2[2] -= 1;
    }
    Real_scalar_field phi2(fftwh.padded_shape_real().vector(),
                           phi_hat2.get_physical_size(),
                           phi_hat2.get_physical_offset(),
                           fftwh.guard_lower(), fftwh.guard_upper());
    //  phi2.get_points().set_storage_size(fftwh.local_size());
    timer("misc");
    fftwh.inv_transform(phi_hat2, phi2);
    timer("invfft phi");
    double norm = 1.0 / (num_points2[0] * num_points2[1] * num_points2[2]);
    phi2.get_points().scale(norm);
    timer("calc norm");
    return phi2;
}

Real_scalar_field
get_phi(Real_scalar_field &rho, Real_scalar_field &phi2, Fftw_helper &fftwh)
{
    // step 6
    Real_scalar_field phi(rho.get_points().get_shape(), rho.get_physical_size(),
                          rho.get_physical_offset(),
                          fftwh.guard_lower(),
                          std::min(fftwh.guard_upper(),
                                   rho.get_points().get_shape()[0]));
    Int3 shape(phi.get_points().get_shape());
    Int3 point, p0;
    timer("misc");
    int i_max = std::min(fftwh.upper(), shape[0]);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (int i = fftwh.lower(); i < i_max; ++i) {
        point[0] = i;
        for (int j = 0; j < shape[1]; ++j) {
            point[1] = j;
            for (int k = 0; k < shape[2]-1; ++k) {
             point[2] = k;
                   phi.get_points().set(point, phi2.get_points().get(point));                     		
            }
	    point[2]=shape[2]-1;
	    p0[0]=i;
	    p0[1]=j;
	    p0[2]=0;
	    phi.get_points().set(point, phi2.get_points().get(p0));
        }
    }


    timer("calc phi");
    return phi;
}



Real_scalar_field
solver_fftw_open(Real_scalar_field &rho, Fftw_helper &fftwh, bool z_periodic,
    bool use_guards)
{
    // The plan: Solve del^2 phi = rho by:
    //  1) convert rho to rho2, where 2 suffix indicates
    //     that we are using Hockney's doubled grid. rho2 = zero when
    //     index outside of boundaries of rho (and so on for other vars).
    //  2) get FFT(rho2) = rho_hat2
    //  3) get (Green function) G_hat2 (intrinsically complex)
    //  4) calculate phi_hat2 = rho_hat2 * G_hat2
    //  5) calculate phi2 = invFFT(phi_hat2)
    //  6) extract phi from phi2


    


    reset_timer();
    //~ Fftw_helper fftwh(rho);
    //~ timer("get plans");
    gather_rho(rho, fftwh.upper());
    timer("gather rho");
    Complex_scalar_field rho_hat2 = get_rho_hat2(rho, fftwh);
    Complex_scalar_field G_hat2 = get_G_hat2(rho, z_periodic, fftwh);
    Complex_scalar_field phi_hat2 = get_phi_hat2(rho, rho_hat2, G_hat2, fftwh);
    timer("misc");
    Real_scalar_field phi2 = get_phi2(rho, phi_hat2, fftwh, z_periodic);
    timer("misc");
    Real_scalar_field phi = get_phi(rho, phi2, fftwh);
 

    timer("misc");
    if (use_guards) {
        fill_guards(phi, fftwh);
    }
    timer("fill guards");
    return phi;
}

