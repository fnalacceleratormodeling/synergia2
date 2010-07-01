#include <mpi.h>
#include "cylindrical.h"
#include "math_constants.h"
#include <fftw3.h>
#include <cmath>
#include "basic_toolkit/PhysicsConstants.h"
#include "mytimer.h"
#include <fstream>


void
get_cylindrical_coords(Macro_bunch_store &mbs, Array_2d<double> &coords)
{
    
  //  boost::multi_array<double, 2> coords_temp(boost::extents[6][mbs.local_num]);
  //  boost::multi_array_types::index_gen index;
    for (int n = 0; n < mbs.local_num; ++n) {
        double x = mbs.local_particles(0,n);
        double y = mbs.local_particles(2,n);
        double r = sqrt(x*x + y*y);
        double theta;
      //  if (r == 0.0) {
        if (r < 1.e-20) {
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
    
//     boost::multi_array_ref<double, 2> coords_temp_ref(coords.get_data_ptr(),boost::extents[6][mbs.local_num]);
//     boost::multi_array<double, 2> coords_temp(coords_temp_ref);
//     for (int n = 0; n < mbs.local_num; ++n) {
//         std::cout<<"coords="<<coords(1,n)
//                 <<" coords_temp_ref="<<coords_temp_ref[1][n]
//                 <<" coords_temp="<<coords_temp[1][n]
//                 <<std::endl;
//    }    
    
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
    std::vector<int> shape(fdomain.get_grid_shape());
//    std::vector<double> cell_size(fdomain.get_cell_size());


    rho.set_all(0.0);

    double cell_r=fdomain.get_cell_size()[0];
    double cell_phi=fdomain.get_cell_size()[1];
    double cell_z=fdomain.get_cell_size()[2];
    bool z_periodic=fdomain.get_periodic_z();
    double total_charge_per_cell_vol=0.;
    for (int n = 0; n < mbs.local_num; ++n) {
        double r = coords(0,n);
        double theta = coords(1,n);
        double z = coords(2,n);
        fdomain.get_leftmost_indices_offsets(r,theta,z,indices,offsets);
        double weight0;


       if (indices[0]>=0) {
          weight0= 1.0 / ((indices[0]+1.)*cell_r*cell_r*cell_phi*cell_z);

	   for (int i = 0; i < 2; ++i) {
              for (int j = 0; j < 2; ++j) {
                 for (int k = 0; k < 2; ++k) {
                    double weight = weight0 * (1 - i - (1 - 2 * i) * offsets[0]) *
                                    (1 - j - (1 - 2 * j) * offsets[1]) *
                                    (1 - k - (1 - 2 * k) * offsets[2]);
                    int ir=indices[0] + i;
                    int iphi=indices[1] + j;
                    int iz=indices[2] + k;
	            if (iphi==shape[1]) iphi=0;
                    if ((iz==shape[2]) && (z_periodic)) iz=0;

                 if (rho.bounds_check(ir,iphi,iz)) {
                            rho(ir,iphi,iz) += weight;
                            total_charge_per_cell_vol += weight*((indices[0]+1.)*cell_r*cell_r*cell_phi*cell_z);
                    }

                  }
              }
          }
      }
      else if (indices[0]==-1){
          //std::cout<<" indices="<<indices[0]<<"  "<<indices[1]<<"  "<<indices[2]<<"  "<<std::endl;
          weight0= 1.0/(pi*(0.5*cell_r)*(0.5*cell_r)*cell_z);
          int ir=0;
          if (r >=1.e-20) {
             for (int j = 0; j < 2; ++j) {
                 for (int k = 0; k < 2; ++k) {
                 double weight =weight0 *(1 - j - (1 - 2 * j) * offsets[1]) *
                                    (1 - k - (1 - 2 * k) * offsets[2]);
                    int iz=indices[2] + k;
                    int iphi=indices[1] + j;
                    if ((iz==shape[2]) && (z_periodic)) iz=0;
                    if (iphi==shape[1]) iphi=0;
                    if (rho.bounds_check(ir,iphi,iz))  {
                        rho(ir,iphi,iz) += weight;
                      total_charge_per_cell_vol += weight*pi*(0.5*cell_r)*(0.5*cell_r)*cell_z;
                    }
                }
            }
         } // AM: to be corrected, put r=0 distribution even on the phi grid, however small effect expected

      }
      else {std::cout<<"error on the grid, indices[0]="<<indices[0]<<" r="<<r<<std::endl;}

    }

    //  std::cout<<"total charge= "<<total_charge_per_cell_vol<<" total number="<<mbs.local_num<<std::endl;


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




     alpha(0) = diag(0);
     z(0) = rhs(0);


     for (i = 1; i < N; ++i) {
           const std::complex<double> t = belowdiag(i)/alpha(i-1);
           alpha(i) = diag(i) - abovediag(i - 1)*t;    //*belowdiag(i)/alpha(i-1);

           if (alpha(i) == 0.0) {
               throw
                std::runtime_error("solve_tridiag_nonsym: zero pivot");
            }

           z(i) = rhs(i) - z(i-1) *t;// *belowdiag(i)/alpha(i-1);

     }


	/* Now back substitute. */
       x(N - 1) = z(N - 1)/alpha(N - 1);
      if (N >= 2) {
         for (i = N - 2, j = 0; j <= N - 2; ++j, --i) {
          x(i) = (z(i) - abovediag(i) * x(i + 1))/alpha(i);
          }
      }


}

void solve_tridiag_nonsymb 
(const boost::multi_array<std::complex<double>, 1> &diag,
         const boost::multi_array<std::complex<double>, 1> &abovediag,
         const boost::multi_array<std::complex<double>, 1> &belowdiag,
         const boost::multi_array<std::complex<double>, 1>::array_view<1>::type &rhs,
         boost::multi_array<std::complex<double>, 1>::array_view<1>::type &x) 

{

     
//     for ( boost::multi_array<std::complex<double>, 1>::iterator iter=x.begin(); iter != x.end(); ++iter){
//         *iter=std::complex<double>(0.,0.);
//     }
    
    for (  int i=0; i<x.size();  ++i){
        x[i]=std::complex<double>(0.,0.);
    }
    
    
    int N = diag.shape()[0];
     
     boost::multi_array<std::complex<double>, 1> alpha(boost::extents[N]);
     boost::multi_array<std::complex<double>, 1>  z(boost::extents[N]);
     size_t i, j;




    alpha[0] = diag[0];
    z[0] = rhs[0];


    for (i = 1; i < N; ++i) {
        const std::complex<double> t = belowdiag[i]/alpha[i-1];
        alpha[i] = diag[i] - abovediag[i - 1]*t;    //*belowdiag(i)/alpha(i-1);

        if (alpha[i] == 0.0) {
            throw
                    std::runtime_error("solve_tridiag_nonsym: zero pivot");
        }

        z[i] = rhs[i] - z[i-1] *t;// *belowdiag(i)/alpha(i-1);

    }


    /* Now back substitute. */
    x[N - 1] = z[N - 1]/alpha[N - 1];
    if (N >= 2) {
        for (i = N - 2, j = 0; j <= N - 2; ++j, --i) {
            x[i] = (z[i] - abovediag[i] * x[i + 1])/alpha[i];
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
decompose_1d(int length, int processors, std::vector<int> &offsets,
		std::vector<int> &counts)
{
	int min_counts = length/processors;
//	int remainder = fmod(length,processors);
    int remainder = length%processors;
	int offset = 0;
	for(int i=0; i<processors; ++i) {
		int count = min_counts;
		if (i >= (processors - remainder)) {
			count += 1;
		}
		offsets.at(i) = offset;
		counts.at(i) = count;
		offset += count;
	}
}


void
solve_cylindrical_finite_periodic(const Cylindrical_field_domain &fdomain,
                                          Array_3d<double > &local_rho, Array_3d<double> &phi)
{

    
    
    timer("start");
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    typedef boost::multi_array<double, 3> array3D_type;
    typedef boost::multi_array_ref<double, 3> array3D_ref_type;
    typedef boost::multi_array<std::complex<double>, 3> array3Dc_type;
    typedef boost::multi_array_ref<std::complex<double>, 3> array3Dc_ref_type;
    typedef boost::multi_array<std::complex<double>, 1> array1Dc_type;
    typedef boost::multi_array_ref<std::complex<double>, 1> array1Dc_ref_type;
    typedef boost::multi_array_types::index_range range;
    
  
      
    array3D_type rho(boost::extents[local_rho.get_shape()[0]][local_rho.get_shape()[1]][local_rho.get_shape()[2]]);
    MPI_Allreduce(reinterpret_cast<void*>(local_rho.get_data_ptr()),
    reinterpret_cast<void*>(rho.data()),
    local_rho.get_size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    timer("gather rho");
    
    std::vector<int> shape = fdomain.get_grid_shape();
    double z_length = fdomain.get_length();
    // the shape of the FFT'd array (shape_lm) is halved in the third
    // dimension because of the peculiar (but efficient) way FFTW does
    // real-to-complex transforms.
    boost::array<array3Dc_type::index, 3> shape_lm = {{shape[0],shape[1],shape[2]/2+1}};   
    array3Dc_type rho_lm(shape_lm);
    
    std::vector<int> offsets(size), counts(size);
    decompose_1d(shape[0],size,offsets,counts);

  
    array3Dc_type sliced_rho_lm(boost::extents[counts[rank]][shape_lm[1]][shape_lm[2]]);
    
  
    if (counts[rank]>0) {
        array3D_type::element* origin_data=rho.data()+offsets[rank]*rho.strides()[0];
       // array3Dc_type::element* origin_fft_data= sliced_rho_lm.data();      
         fftw_plan plan = fftw_plan_many_dft_r2c(2,
                &shape[1], counts[rank],
                origin_data, // sliced_rho_ref.data(),                   
                NULL, 1, shape[1]*shape[2],
                reinterpret_cast<double (*)[2]>(sliced_rho_lm.data()),
                NULL, 1, shape_lm[1]*shape_lm[2],
                FFTW_ESTIMATE);
        timer("plan");
        fftw_execute(plan);
        timer("fft");

    }
    
    std::vector<int> receive_offsets(size), receive_counts(size);
    for (int i=0; i< size; ++i) {
        receive_counts.at(i) = counts.at(i)*shape_lm[1]*shape_lm[2];
        receive_offsets.at(i) = offsets.at(i)*shape_lm[1]*shape_lm[2];
    }

//     MPI_Allgatherv(reinterpret_cast<void*>(sliced_rho_lm.data()),
//                    counts[rank]*shape_lm[1]*shape_lm[2],
//                    MPI_DOUBLE_COMPLEX,
//                    reinterpret_cast<void*>(rho_lm.data()),
//                   &receive_counts[0],
//                   &receive_offsets[0],
//                    MPI_DOUBLE_COMPLEX,
//                    MPI_COMM_WORLD);
//  AM: I suspect a MPI bug in MPI_Allgatherv....therefore I replace  MPI_Allgatherv with a loop of MPI_Gatherv which seems to work fine   
    for (int rnk=0;rnk<size;++rnk){
        MPI_Gatherv(reinterpret_cast<void*>(sliced_rho_lm.data()),
                       counts[rank]*shape_lm[1]*shape_lm[2],
                       MPI_DOUBLE_COMPLEX,
                       reinterpret_cast<void*>(rho_lm.data()),
                                               &receive_counts[0],
                                               &receive_offsets[0],
                                               MPI_DOUBLE_COMPLEX,
                                               rnk,
                                               MPI_COMM_WORLD);       
    } 
    timer("gather");
    std::vector<int> offsets2(size), counts2(size);
    decompose_1d(shape_lm[1],size,offsets2,counts2);
    
    array3D_type::index_gen indices;  
    array1Dc_type diag(boost::extents[shape_lm[0]]);
    array1Dc_type  above_diag(boost::extents[shape_lm[0]]);
    array1Dc_type below_diag(boost::extents[shape_lm[0]]);
    array3Dc_type phi_lm(shape_lm),phi_lm_local(shape_lm);
    
    double deltar = fdomain.get_radius()/(fdomain.get_grid_shape()[0]+0.5);
    double r;
    for(int l=offsets2[rank]; l<(offsets2[rank]+counts2[rank]); ++l) { 
        int wavenumber_l = (l+shape_lm[1]/2) % shape_lm[1] -
                shape_lm[1]/2;
        for(int m=0; m<shape_lm[2]; ++m) {
            int wavenumber_m = (m+shape_lm[2]/2) % shape_lm[2] -
                    shape_lm[2]/2;
            for(int i=0; i<shape_lm[0]; ++i) {            
                r = (i+0.5)*deltar;
                diag[i] = -2.0/(deltar*deltar);
                diag[i] += - wavenumber_l*wavenumber_l/(r*r)
                        -(2*pi*wavenumber_m/z_length)*(2*pi*wavenumber_m/z_length);
                if (i<(shape_lm[0]-1)) {
                    above_diag[i] = 1.0/(deltar*deltar) + 1.0/(2*deltar*r);
                } else {above_diag[i]=std::complex<double>(0.,0.);}
                if (i>0) {
                    below_diag[i] = 1.0/(deltar*deltar) - 1.0/(2*deltar*r);
                }else {below_diag[i]=std::complex<double>(0.,0.);}
            }
            array1Dc_type::array_view<1>::type rhs=rho_lm[indices[range()][l][m]];
            array1Dc_type::array_view<1>::type x=phi_lm_local[indices[range()][l][m]];
            solve_tridiag_nonsymb(diag,above_diag,below_diag,rhs,x);
        }
    }
    timer("tridiag");
    std::vector<int> receive_offsets2(size), receive_counts2(size);
    for (int i=0; i< size; ++i) {
        receive_counts2.at(i) = counts2.at(i)*shape_lm[2];
        receive_offsets2.at(i) = offsets2.at(i)*shape_lm[2];
    }
    
    for(int i=0; i<shape_lm[0]; ++i) {
      //  if  (counts2[rank]>0){
            array3Dc_type::element* send_origin=
                phi_lm_local.origin()+i*phi_lm_local.strides()[0]+offsets2[rank]*phi_lm_local.strides()[1];
            array3Dc_type::element* receive_origin=
                 phi_lm.origin()+i*phi_lm.strides()[0];
       
//        MPI_Allgatherv(reinterpret_cast<void*>(send_origin),                                           
//                       counts2[rank]*shape_lm[2],
//                        MPI_DOUBLE_COMPLEX,
//                       reinterpret_cast<void*>(receive_origin),
//                         &receive_counts2[0],
//                         &receive_offsets2[0],
//                          MPI_DOUBLE_COMPLEX,
//                          MPI_COMM_WORLD); 
//  AM: I suspect a MPI bug in MPI_Allgatherv....therefore I replace  MPI_Allgatherv with a loop of MPI_Gatherv which seems to work fine             
            for (int rnk=0;rnk<size;++rnk){
                MPI_Gatherv(reinterpret_cast<void*>(send_origin),                                           
                         counts2[rank]*shape_lm[2],
                         MPI_DOUBLE_COMPLEX,
                         reinterpret_cast<void*>(receive_origin),
                         &receive_counts2[0],
                         &receive_offsets2[0],
                         MPI_DOUBLE_COMPLEX,
                         rnk,
                         MPI_COMM_WORLD);    
            } 
       // }  
    }
    timer("gather2");
  
    // FFTW transforms are not normalized. We need to apply the normalization
    // manually.      
    array3D_ref_type phi_ref(phi.get_data_ptr(),
                                           boost::extents[shape[0]][shape[1]][shape[2]]);

    
     if (counts[rank]>0) {
        fftw_plan plani = fftw_plan_many_dft_c2r(2,
                &shape[1], counts[rank],
                reinterpret_cast<double (*)[2]>(phi_lm.data()+offsets[rank]*phi_lm.strides()[0]),
                        NULL, 1, shape_lm[1]*shape_lm[2],
                        phi_ref.data()+offsets[rank]*phi_ref.strides()[0],                        
                                              NULL, 1, shape[1]*shape[2],
                                              FFTW_ESTIMATE);
        timer("invplan");
        fftw_execute(plani);
        timer("ifft");
    }

    /*
    phi.slice(
              Range(offsets[rank],offsets[rank]+counts[rank]-1),
                    Range(),Range(),false).scale(1.0/(shape[1]*shape[2]));*/
        
       
    for (int i=offsets[rank]; i<offsets[rank]+counts[rank]; ++i){
        for (int j=0; j<shape[1];++j) {
            for (int k=0; k<shape[2];++k) {
                phi_ref[i][j][k] *=  1.0/(shape[1]*shape[2]); 
            }
        }
    }
        
   
  
    // at end of routine, phi is the size of the global array, but it only contains
    // the locally calculated pieces.
}




void
fill_guards_cylindrical(const Cylindrical_field_domain &fdomain,
		Array_3d<double> &phi)
{
	std::vector<int> shape = fdomain.get_grid_shape();
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<int> offsets(size), counts(size);
	decompose_1d(shape[0],size,offsets,counts);

	void *send_buffer, *recv_buffer;
    
	MPI_Status status;
 	int message_size = shape[1] * shape[2];
	int lower = offsets[rank];
	int upper = offsets[rank] + counts[rank]-1;
    
	// send to the right
    typedef boost::multi_array_ref<double, 3> array3D_ref_type;
    array3D_ref_type phi_ref(phi.get_data_ptr(),
                             boost::extents[shape[0]][shape[1]][shape[2]]);
    
    if ((upper < shape[0]-1) &&  (upper >= 0)) {
        send_buffer=reinterpret_cast<void*>(phi_ref.data()+upper*phi_ref.strides()[0]);
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, rank + 1, rank,
        		MPI_COMM_WORLD);
    }
    if (lower > 0) {
        recv_buffer=reinterpret_cast<void*>(phi_ref.data()+(lower-1)*phi_ref.strides()[0]);
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, rank - 1, rank - 1,
        		MPI_COMM_WORLD, &status);
    }

	// send to the left
	if (lower > 0) {
        send_buffer=reinterpret_cast<void*>(phi_ref.data()+lower*phi_ref.strides()[0]);
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, rank - 1, rank,
        		MPI_COMM_WORLD);
    }
    if ((upper < shape[0]-1) &&  (upper >= 0)){
        recv_buffer=reinterpret_cast<void*>(phi_ref.data()+(upper+1)*phi_ref.strides()[0]);        
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, rank + 1, rank + 1,
        		MPI_COMM_WORLD, &status);
    }
}

void
calculate_E_cylindrical(const Cylindrical_field_domain &fdomain,
                          Array_3d<double> &phi,
                          Array_3d<double> &Ex,
                          Array_3d<double> &Ey,
                          Array_3d<double> &Ez)
{ 
    
   
	timer("misc calc_E");
	fill_guards_cylindrical(fdomain,phi);
	timer("fill guards");
    std::vector<int> shape = phi.get_shape();
    
    typedef boost::multi_array_ref<double, 3> array3D_ref_type;
    array3D_ref_type Ex_ref(Ex.get_data_ptr(),
                            boost::extents[shape[0]][shape[1]][shape[2]]);
    
    array3D_ref_type Ey_ref(Ey.get_data_ptr(),
                            boost::extents[shape[0]][shape[1]][shape[2]]);
    
    
    array3D_ref_type Ez_ref(Ez.get_data_ptr(),
                            boost::extents[shape[0]][shape[1]][shape[2]]);
    
    array3D_ref_type phi_ref(phi.get_data_ptr(),
                            boost::extents[shape[0]][shape[1]][shape[2]]);
    
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::vector<int> offsets(size), counts(size);
	std::vector<int> receive_offsets(size), receive_counts(size);
	decompose_1d(shape[0],size,offsets,counts);
	for (int i=0; i< size; ++i) {
		receive_counts.at(i) = counts.at(i)*shape[1]*shape[2];
		receive_offsets.at(i) = offsets.at(i)*shape[1]*shape[2];
	}
  
    double theta_step = 2.0*pi/shape[1];
    double z_step = fdomain.get_length()/shape[2];
    double ordinary_r_step = fdomain.get_radius()/(shape[0]+0.5);
    for(int i_r = offsets[rank]; i_r<offsets[rank]+counts[rank]; ++i_r) { 
        int r_left = i_r-1;
        int r_right = i_r+1;
        double r_step = 2.*ordinary_r_step;
        double r = (i_r+0.5)*ordinary_r_step;
         if(i_r == 0) {
             r_left = 0;
             r_step *= 0.5;
         }

        for(int i_theta = 0; i_theta<shape[1]; ++i_theta) {
            int theta_left = i_theta-1;
            int theta_right = i_theta+1;
            double theta = i_theta*theta_step;
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
                     z_left = shape[2]-1;
                 }
                 if (z_right == shape[2]) {
                     z_right = 0;
                 }

// 	 if (!phi.bounds_check(r_right,i_theta,i_z)) { std::cout<<"out of bounds phi1"<<std::endl;
//          std::cout<<"r_right="<<r_right<<" i_theta="<<i_theta<<"  i_z="<<i_z<<std::endl;; abort();}
//          if (!phi.bounds_check(r_left,i_theta,i_z)) { std::cout<<"out of bounds phi2"<<std::endl;  abort();}
//          if (!phi.bounds_check(i_r,theta_right,i_z)) { std::cout<<"out of bounds phi3"<<std::endl;  abort();}
//          if (!phi.bounds_check(i_r,theta_left,i_z)) { std::cout<<"out of bounds phi4"<<std::endl;  abort();}
//          if (!phi.bounds_check(i_r,i_theta,z_right)) { std::cout<<"out of bounds phi5"<<std::endl;
//                   std::cout<<"i_r="<<i_r<<" i_theta="<<i_theta<<"  z_right="<<z_right<<std::endl;
//                                   abort();}
//          if (!phi.bounds_check(i_r,i_theta,z_left)) { std::cout<<"out of bounds phi6"<<std::endl;  abort();}



	       double dphi_dr;
              // if(i_r == shape[0]-1) {dphi_dr=- phi(r_left,i_theta,i_z)/(r_step); // because phi(shape[0],i_theta,i_z)=0
               if(i_r == shape[0]-1) {dphi_dr=- phi_ref[r_left][i_theta][i_z]/(r_step); // because phi(shape[0],i_theta,i_z)=0   
                }
                else   {dphi_dr = (phi_ref[r_right][i_theta][i_z] - phi_ref[r_left][i_theta][i_z])/r_step;
                }


                double dphi_dtheta = (phi_ref[i_r][theta_right][i_z] -
                            phi_ref[i_r][theta_left][i_z])/(2.0*theta_step);
                double dphi_dz = (phi_ref[i_r][i_theta][z_right] -
                            phi_ref[i_r][i_theta][z_left])/(2.0*z_step);

/*
            if ((!Ex.bounds_check(i_r,i_theta,i_z)) || (!Ey.bounds_check(i_r,i_theta,i_z)) || (!Ez.bounds_check(i_r,i_theta,i_z)) ) { std::cout<<"out of bounds E"<<std::endl;
	    abort();}*/

                Ex_ref[i_r][i_theta][i_z] = cos(theta)*dphi_dr - sin(theta)*dphi_dtheta/r;
                Ey_ref[i_r][i_theta][i_z] = sin(theta)*dphi_dr + cos(theta)*dphi_dtheta/r;
                Ez_ref[i_r][i_theta][i_z] = dphi_dz;

            }
        }
    }
    timer("calcE");
     if (size > 1) {

//         int place= Ex_ref.strides()[0]*offsets[rank];   
//         MPI_Allgatherv(reinterpret_cast<void*>(Ex_ref.data()+place),
//                 receive_counts[rank],
//                 MPI_DOUBLE,
//                 reinterpret_cast<void*>(Ex_ref.data()),
//                 &receive_counts[0],
//                 &receive_offsets[0],
//                 MPI_DOUBLE,
//                 MPI_COMM_WORLD); 
//        
//         MPI_Allgatherv(reinterpret_cast<void*>(Ey_ref.data()+place),
//                 receive_counts[rank],
//                 MPI_DOUBLE,
//                 reinterpret_cast<void*>(Ey_ref.data()),
//                 &receive_counts[0],
//                 &receive_offsets[0],
//                 MPI_DOUBLE,
//                 MPI_COMM_WORLD);
//        
//         MPI_Allgatherv(reinterpret_cast<void*>(Ez_ref.data()+place),
//                 receive_counts[rank],
//                 MPI_DOUBLE,
//                 reinterpret_cast<void*>(Ez_ref.data()),
//                 &receive_counts[0],
//                 &receive_offsets[0],
//                 MPI_DOUBLE,
//                 MPI_COMM_WORLD);  
//  AM: I suspect a MPI bug in MPI_Allgatherv....therefore I replace  MPI_Allgatherv with a loop of MPI_Gatherv which seems to work fine             
        for (int rnk=0;rnk<size;++rnk){  
            int place= Ex_ref.strides()[0]*offsets[rank];   
            MPI_Gatherv(reinterpret_cast<void*>(Ex_ref.data()+place),
                receive_counts[rank],
                MPI_DOUBLE,
                reinterpret_cast<void*>(Ex_ref.data()),
                &receive_counts[0],
                &receive_offsets[0],
                MPI_DOUBLE,
                rnk,
                MPI_COMM_WORLD); 
       
            MPI_Gatherv(reinterpret_cast<void*>(Ey_ref.data()+place),
                receive_counts[rank],
                MPI_DOUBLE,
                reinterpret_cast<void*>(Ey_ref.data()),
                &receive_counts[0],
                &receive_offsets[0],
                MPI_DOUBLE,
                rnk,
                MPI_COMM_WORLD);
       
            MPI_Gatherv(reinterpret_cast<void*>(Ez_ref.data()+place),
                receive_counts[rank],
                MPI_DOUBLE,
                reinterpret_cast<void*>(Ez_ref.data()),
                &receive_counts[0],
                &receive_offsets[0],
                MPI_DOUBLE,
                rnk,
                MPI_COMM_WORLD);    

        }
        
     }
    
    
    
    timer("gatherE");
 // std::cout<< " end E_cylindrical"<<std::endl;
}



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

// inline
// double
// interpolate_3d_cyl_old(double x1, double x2, double x3,
//                const Cylindrical_field_domain &fdomain,
//                const Array_3d<double> &points)
// {
//     // Interpolate between grid points. There is no unique scheme to do this
//     // in 3D, so we choose to use trilinear interpolation.
//     std::vector<int> c(3); // c for corner
//     std::vector<double> f(3); // f for fractional difference
//     fdomain.get_leftmost_indices_offsets(x1,x2,x3,c,f);
//     double val = ((1.0 - f[0]) * (1.0 - f[1]) * (1.0 - f[2]) * points.at(c[0],c[1],c[2]) +
//             f[0] * (1.0 - f[1]) * (1.0 - f[2]) * points.at(c[0] + 1, c[1], c[2]) +
//             (1.0 - f[0]) * f[1] * (1.0 - f[2]) * points.at(c[0], c[1] + 1, c[2]) +
//             (1.0 - f[0]) * (1.0 - f[1]) * f[2] * points.at(c[0], c[1], c[2] + 1) +
//             f[0] * f[1] * (1.0 - f[2]) * points.at(c[0] + 1, c[1] + 1, c[2]) +
//             f[0] * (1.0 - f[1]) * f[2] * points.at(c[0] + 1, c[1], c[2] + 1) +
//             (1.0 - f[0]) * f[1] * f[2] * points.at(c[0], c[1] + 1, c[2] + 1) +
//             f[0] * f[1] * f[2] * points.at(c[0] + 1, c[1] + 1, c[2] + 1));
//     return val;
// }

inline
double
interpolate_3d_cyl(double x1, double x2, double x3,
               const Cylindrical_field_domain &fdomain,
               const Array_3d<double> &points)
{

    std::vector<int> shape=fdomain.get_grid_shape();
    std::vector<int> c(3); // c for corner
    std::vector<double> f(3); // f for fractional difference
    fdomain.get_leftmost_indices_offsets(x1,x2,x3,c,f);
    int c0=c[0];
    int c1=c[1];
    int c2=c[2];
    int c0p=c0+1;
    int c1p=c1+1;
    int c2p=c2+1;


    if(c1p==shape[1]) c1p=0;
    if(c2p==shape[2]) c2p=0;

//     std::cout<<"shape="<<shape[0]<<"  "<<shape[1]<<"  "<<shape[2]<<std::endl;
//     std::cout<<"c0= "<<c0<<"  c1= "<<c1<<"  c2= "<<c2<<std::endl;
//     std::cout<<"c0p= "<<c0p<<"  c1p= "<<c1p<<"  c2p= "<<c2p<<std::endl;

    double val;

    if ((c0>=0) && (c0<shape[0])) {

// 	    if ((!points.bounds_check(c0,c1,c2)) || (!points.bounds_check(c0,c1p,c2)) || (!points.bounds_check(c0,c1,c2p)) || (!points.bounds_check(c0,c1p,c2p))) { std::cout<<"out of bounds 1"<<std::endl;
// 	    std::cout<<"c0="<<c0<<"  c1="<<c1<<"  c2="<<c2<<"  c1p="<<c1p<<" c2p="<<c2p<<std::endl;
//             std::cout<<std::setprecision(18)<<"x1="<<x1<<"  x2="<<x2<<"  x3="<<x3<<"  0.5 length="<<0.5*fdomain.get_length()<<std::endl;
//             std::cout<<"shape ="<<shape[0]<<" "<<shape[1]<<" "<<shape[2]<<std::endl;
// 	    return 0.;}
            val = ((1.0 - f[0]) * (1.0 - f[1]) * (1.0 - f[2]) * points.at(c0,c1,c2) +
            (1.0 - f[0]) * f[1] * (1.0 - f[2]) * points.at(c0, c1p, c2) +
            (1.0 - f[0]) * (1.0 - f[1]) * f[2] * points.at(c0, c1, c2p) +
            (1.0 - f[0]) * f[1] * f[2] * points.at(c0, c1p, c2p));


          if (c0p<shape[0]){
//                               if ((!points.bounds_check(c0p,c1,c2)) || (!points.bounds_check(c0p,c1p,c2)) || (!points.bounds_check(c0p,c1,c2p)) || (!points.bounds_check(c0p,c1p,c2p))) { std::cout<<"out of bounds 2"<<std::endl;
//                         std::cout<<"c0p="<<c0p<<"  c1="<<c1<<"  c2="<<c2<<"  c1p="<<c1p<<" c2p="<<c2p<<std::endl;
// 			std::cout<<std::setprecision(18)<<"x1="<<x1<<"  x2="<<x2<<"  x3="<<x3<<"  0.5 length="<<0.5*fdomain.get_length()<<std::endl;
//                         std::cout<<"shape ="<<shape[0]<<" "<<shape[1]<<" "<<shape[2]<<std::endl;
//
//
// 	                       return 0.;}

	                    val +=   (f[0] * (1.0 - f[1]) * (1.0 - f[2]) * points.at(c0p, c1, c2) +
                            f[0] * f[1] * (1.0 - f[2]) * points.at(c0p, c1p, c2) +
                            f[0] * (1.0 - f[1]) * f[2] * points.at(c0p, c1, c2p) +
                            f[0] * f[1] * f[2] * points.at(c0p, c1p, c2p));}
     }
     else if (c0==-1) {

// 		if ((!points.bounds_check(0,c1,c2)) || (!points.bounds_check(0,c1p,c2)) || (!points.bounds_check(0,c1,c2p)) || (!points.bounds_check(0,c1p,c2p))) { std::cout<<"out of bounds 3"<<std::endl;
//              std::cout<<"  c1="<<c1<<"  c2="<<c2<<"  c1p="<<c1p<<" c2p="<<c2p<<std::endl;
//             std::cout<<std::setprecision(18)<<"x1="<<x1<<"  x2="<<x2<<"  x3="<<x3<<"  0.5 length="<<0.5*fdomain.get_length()<<std::endl;
//             std::cout<<"shape ="<<shape[0]<<" "<<shape[1]<<" "<<shape[2]<<std::endl;
//
//
//
// 	                       return 0.;}

            val = ((1.0 - f[1]) * (1.0 - f[2]) * points.at(0,c1,c2) +
            f[1] * (1.0 - f[2]) * points.at(0, c1p, c2) +
            (1.0 - f[1]) * f[2] * points.at(0, c1, c2p) +
            f[1] * f[2] * points.at(0, c1p, c2p));
     }
   //  else{ std::cout<<" indices[0] error, in interoplation, indices[0]="<<c[0]<<std::endl;}

    return val;
}

void
full_kick_cylindrical(const Cylindrical_field_domain &fdomain,
                      Array_3d<double> &phi, double tau,
                      Macro_bunch_store &mbs, Array_2d<double> &coords)
{
    std::vector<int> shape = fdomain.get_grid_shape();
    Array_3d<double> Ex(shape[0],shape[1],shape[2]);
    Array_3d<double> Ey(shape[0],shape[1],shape[2]);
    Array_3d<double> Ez(shape[0],shape[1],shape[2]);
    Ex.set_all(0.);
    Ey.set_all(0.);
    Ez.set_all(0.);
    calculate_E_cylindrical(fdomain,phi,Ex,Ey,Ez);



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
    factor *= 1.0/(beta * c*mbs.total_num); // the  arc length tau=beta*c* (Delta t), so (Delta t)= tau/(beta*c)
    factor *= 1.0/gamma;    // relativistic factor,
    factor *=mbs.units(1); // the kikcing force should be muliplied  by the unit of p, this is a factor of 1/mass

    double xyfactor =tau*factor;
    double zfactor =tau*factor*gamma;

    

    for (int n = 0; n < mbs.local_num; ++n) {
        double r = coords(0,n);
        double theta = coords(1,n);
        double z = coords(2,n);


        mbs.local_particles(1,n) += xyfactor*
                interpolate_3d_cyl(r,theta,z,fdomain,Ex);

        mbs.local_particles(3,n) += xyfactor*
                interpolate_3d_cyl(r,theta,z,fdomain,Ey);

        mbs.local_particles(5,n) += zfactor*
                interpolate_3d_cyl(r,theta,z,fdomain,Ez);

    }

}
