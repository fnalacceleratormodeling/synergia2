#include "impedance.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/period.h"
#include "synergia/utils/simple_timer.h"

// Impedance::Impedance(double const & orbit_length, int const & zgrid):
//      Collective_operator("impedance"), orbit_length(orbit_length), z_grid(zgrid), wake_file("no_file")
// {  
//      wake_factor=-4.*mconstants::pi*pconstants::rp/pconstants::c;
// }

Impedance::Impedance(std::string const & wake_file, double const & orbit_length, double const & bunchsp, int const & zgrid, std::string const & pipe_symmetry, int const nstored_turns):
     Collective_operator("impedance"), orbit_length(orbit_length), bunch_spacing(bunchsp), z_grid(zgrid), wake_file(wake_file),
     nstored_turns(nstored_turns)
{  
    // wake_factor=-4.*mconstants::pi*pconstants::rp/pconstants::c;
     wake_factor=-4.*mconstants::pi*pconstants::rp;   
     std::cout<<" in impedance wake_factor="<<wake_factor<<std::endl;
     this->pipe_symmetry=pipe_symmetry;
     
     
     
     
     
//  read the wakes from the file wake_file     
// for parallel plates geometry wake file should be written as a four column file such, containing wakes functions such:
//  z[m]        Wz_trx/Z_0/L[1/(m^2*s]       Wz_try/Z_0/L[1/(m^2*s]        Wz_l/Z_0/L[1/(ms)]   
// the line starting with "#" in the file are skipped and can be  used for comments 
     std::ifstream rfile;
     std:: string line;
     rfile.open(wake_file.c_str());
      while (!rfile.eof() && rfile.is_open()) {
          getline(rfile,line);       
          if (line.substr(0,1) != "#" && !line.empty() ){
                double column1, column2, column3, column4;
                char * endword1, * endword2;
                column1=strtod(line.c_str(), &endword1);
                column2=strtod(endword1, &endword2);    
                if (endword1 == endword2 ) { std::cout<<" the wake file should have at least 3 columns, z, W_transverse and W_z"<<std::endl;
                                                abort();               
                                            } 
                // column2 read succesfully                                   
                column3=strtod(endword2, &endword1);
                if (endword1 == endword2  ) { std::cout<<" the wake file should have at least 3 columns, z, W_transverse and W_z"<<std::endl;
                                                abort();               
                                            }                                          
                // column3 read succesfully   
                column4=strtod(endword1, &endword2);                            
                if (endword1 == endword2  )  { // there is no column 4, but can be OK  for circular pipe
                        if (get_pipe_symmetry()=="circular") { 
                            z_coord.push_back(column1);
                            x_wake.push_back(column2);
                            y_wake.push_back(column2);
                            z_wake.push_back(column3);
                        }  
                        else {  std::cout<< " pipe symmetry is "<<get_pipe_symmetry()<<std::endl;
                                std::cout<<" the wake file should have 4 columns, z, W_x, W_y and W_z"<<std::endl;
                                                abort(); 
                        }
                 } 
                 else {
                        z_coord.push_back(column1);
                        x_wake.push_back(column2);
                        y_wake.push_back(column3);
                        z_wake.push_back(column4);
                        
                }                                                          
         
          }    
      }
      rfile.close();
      
      std::cout<<"  wake read from  "<<wake_file<<std::endl;
      std::cout<<"  pipe symmetry  "<<get_pipe_symmetry()<<std::endl;
      std::cout<<"  number of previous turns considered is  "<<get_nstored_turns()<<std::endl;
     // wakes read!
     
     
}

//******************************************************************************

void
calculate_moments_and_partitions(Bunch & bunch, MArray1d &zdensity,  MArray1d &xmom, MArray1d &ymom,  double  z_left, double z_length, MArray1int  &bin_partition)
{

  //  std::cout<<"begining of calculate_moments_and_partitions"<<std::endl;
    Commxx comm(bunch.get_comm());
    int rank(comm.get_rank());
    
  
      int z_num = zdensity.num_elements();
      double h = z_length/(z_num-1.0); // AM , before was h = z_length/z_num
      if (z_length<= 1.e-14 )   throw
                 std::runtime_error("h and z_length too small ");
//     double h = z_length/z_num;  // AM, what if z_periodic??? 


    
     
     MArray1d  local_zdensity(boost::extents[z_num]);
     MArray1d  local_xmom(boost::extents[z_num]);
     MArray1d  local_ymom(boost::extents[z_num]);
     
      
    for (  int i=0; i<z_num;  ++i){
        local_zdensity[i]=0.0;
        local_xmom[i]=0.0;
        local_ymom[i]=0.0;
    }



     for (int part = 0;  part < bunch.get_local_num(); ++part) {
         int bin = static_cast<int>((bunch.get_local_particles()[part][4]-z_left)/h);
         if ((bin < z_num) && (bin >= 0)) {
             local_zdensity[bin] += 1;
             local_xmom[bin] += bunch.get_local_particles()[part][0];
             local_ymom[bin] += bunch.get_local_particles()[part][2];
             bin_partition[part]=bin; //bin_partition(n) is the bin where you find the particle n
         }
         else if ((bin==z_num) && fabs(bunch.get_local_particles()[part][4]-z_length-z_left)<z_length*1.e-14) { 
             local_zdensity[bin-1] += 1; // put the edge particle in the last bin=z_num-1
             bin_partition[part]=z_num-1;
         } 
         else
         {   std::cout << "  z_left  "<<z_left<<"  rank= "<<rank<<std::endl;
          std::cout<<"bunch.get_local_particles()[part][4]="  <<bunch.get_local_particles()[part][4]<<"  rank= "<<rank<<std::endl; 
           std::cout<<"bunch.get_local_particles()[part]0,1,2,3,4,5="  <<bunch.get_local_particles()[part][0]<<
           "  "<<bunch.get_local_particles()[part][1]<<
           "  "<<bunch.get_local_particles()[part][2]<<
           "  "<<bunch.get_local_particles()[part][3]<<
           "  "<<bunch.get_local_particles()[part][4]<<
           "  "<<bunch.get_local_particles()[part][5]<<std::endl;
           
         std::cout<< "NNNNNNNNNN part="<<part<<std::endl; 
         std::cout << "  z_length  "<<z_length<<"  rank= "<<rank<<std::endl;
         std::cout << "(mbs.local_particles(4,n)-z_left)= "<<(bunch.get_local_particles()[part][4]-z_left)<<"  rank= "<<rank<<std::endl;
         std::cout << "bin: " << bin<<"  z_num="<<z_num<< "  h=" << h <<"  rank= "<<rank<<std::endl;                
         std::cout << "bunch.get_local_particles()[part][4]-z_length-z_left= "<<fabs(bunch.get_local_particles()[part][4]-z_length-z_left)<<"  rank= "<<rank<<std::endl;
         throw
                 std::runtime_error("particles out of range in calculate_moments_and_partitions ");
         }

    
    }


    int error = MPI_Allreduce(reinterpret_cast<void*>(local_zdensity.origin()), reinterpret_cast<void*> (zdensity.origin()),z_num, MPI_DOUBLE,
            MPI_SUM, comm.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance zdensity");
    }
    
    
    MPI_Allreduce(reinterpret_cast<void*>(local_xmom.origin()),
                   reinterpret_cast<void*>(xmom.origin()),
                                           z_num, MPI_DOUBLE, MPI_SUM, comm.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance xmom");
    }
  
    
    MPI_Allreduce(reinterpret_cast<void*>(local_ymom.origin()),
                   reinterpret_cast<void*>(ymom.origin()),
                                           z_num, MPI_DOUBLE, MPI_SUM, comm.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Impedance ymom");
    }    

//     //    dbg is set here XXXXXXXXXXXXXX
    int dbg = 0;
  //  std::cout<<"end of moments and partitions"<<std::endl;
    
    for (int k = 0; k < z_num; ++k) {
       // std::cout<<"zdensity[k]="<<zdensity[k]<<std::endl;
        if (zdensity[k] != 0.0) {
            if (dbg) std::cout << "before bin: " << k << " zdensity(k): " << zdensity[k] << " xmom(k): " << xmom[k] << " ymom(k): " << ymom[k] << std::endl;
            xmom[k] /= zdensity[k];
            ymom[k] /= zdensity[k];
            if (dbg) std::cout << "after bin: " << k << " zdensity(k): " << zdensity[k] << " xmom(k): " << xmom[k] << " ymom(k): " << ymom[k] << std::endl;
        } else {
            xmom[k] = 0.0;
            ymom[k] = 0.0;
        }
    }
    
}   
//*************************************************************************
 void 
 get_kicks(int z_grid, double line_length,double  N_factor, double  cell_size_z, std::vector<double> & zcoord, std::vector<double> & xwake, std::vector<double> & ywake, std::vector<double> & zwake, MArray1d & zdensity, MArray1d &  xmom, MArray1d & ymom, double & bunch_sp, std::list<Bunch_means> & stored_bunches, MArray1d & dipole_x, MArray1d &  dipole_y, MArray1d & quad_y, MArray1d & l_monopole)
 {
    
    int zpoints=zcoord.size();
    int registered_turns=stored_bunches.size();
  //  std::cout<<" registred turns= "<<registered_turns<<std::endl; 
 
    for (int i = 0; i < z_grid; ++i){
      // in-bunch impedance  
        for (int j = i+1; j < z_grid; ++j){
            double zji=(j-i)*cell_size_z;
            if (zji>=zcoord[0]) {// at small distance no impedance is considered
                // below it is assumed the wake function is stored using a quadratic grid
                int iz=static_cast<int>(floor(sqrt((zji-zcoord[0])/(zcoord[1]-zcoord[0]))));
                double wake_x(0.), wake_y(0.), wake_z(0.);
                    if (iz+1 <= zpoints) {
                    wake_x=xwake[iz]+(zji-zcoord[iz])*(xwake[iz+1]-xwake[iz])/(zcoord[iz+1]-zcoord[iz]);
                    wake_y=ywake[iz]+(zji-zcoord[iz])*(ywake[iz+1]-ywake[iz])/(zcoord[iz+1]-zcoord[iz]);
                    wake_z=zwake[iz]+(zji-zcoord[iz])*(zwake[iz+1]-zwake[iz])/(zcoord[iz+1]-zcoord[iz]); 
                    }
                dipole_x[i] += zdensity[j]*N_factor*xmom[j]*wake_x; 
                dipole_y[i] += zdensity[j]*N_factor*ymom[j]*wake_y; 
                quad_y[i] += zdensity[j]*N_factor*wake_x;               
                l_monopole[i] += zdensity[j]*N_factor*wake_z; 
            }          
         }
    
    }
    
//  impedance contribution from previous turns    
    double xdipole=0.;
    double ydipole=0.;
    double quadcontrib=0.; 
    double  lmonocontrib=0.;
    
    std::list<Bunch_means>::const_iterator it = stored_bunches.begin();
    for (int iturn=1; iturn<registered_turns; ++iturn){    
            it ++;
            double zji=iturn*line_length; 
            int iz=static_cast<int>(floor(sqrt((zji-zcoord[0])/(zcoord[1]-zcoord[0]))));
            double wake_x(0.), wake_y(0.), wake_z(0.);
            if (iz+1 <= zpoints) {                
               wake_x=xwake[iz]+(zji-zcoord[iz])*(xwake[iz+1]-xwake[iz])/(zcoord[iz+1]-zcoord[iz]);
               wake_y=ywake[iz]+(zji-zcoord[iz])*(ywake[iz+1]-ywake[iz])/(zcoord[iz+1]-zcoord[iz]);
               wake_z=zwake[iz]+(zji-zcoord[iz])*(zwake[iz+1]-zwake[iz])/(zcoord[iz+1]-zcoord[iz]);               
            } 
             xdipole += N_factor*(*it).n_part*(*it).x_mean*wake_x; 
             ydipole += N_factor*(*it).n_part*(*it).y_mean*wake_y; 
             quadcontrib +=N_factor*(*it).n_part*wake_x;
             lmonocontrib += N_factor*(*it).n_part*wake_z; 
        
           //  xdipole += stored_bunchnp(ibunch,iturn)* stored_means(ibunch,iturn,0)*wake_x; 
//             ydipole += stored_bunchnp(ibunch,iturn)* stored_means(ibunch,iturn,1)*wake_y;
//             quadcontrib += stored_bunchnp(ibunch,iturn)*wake_x;
//             lmonocontrib +=  stored_bunchnp(ibunch,iturn)*wake_z;
        
           
    }
     it++;
  //   if (it != stored_bunches.end()) std::cout<<" aaaaaaaaaaaaaaa"<<std::endl;
    
    for (int i = 0; i < z_grid; ++i){
        dipole_x[i] +=xdipole;
        dipole_y[i] +=ydipole;
        quad_y[i] +=quadcontrib;
        l_monopole[i] +=lmonocontrib; 
       }    
 }



void 
impedance_kick(Bunch & bunch, double wake_factor,  MArray1int & bin_partition, MArray1d & dipole_x, MArray1d & dipole_y, MArray1d & quad_y,  MArray1d & l_monopole)
{
   
 for (int part = 0; part < bunch.get_local_num(); ++part) {
        double xkick=0., ykick=0., zkick=0.;
        int bin=bin_partition[part];  // bin_partition(n) is the bin where you find the particle n 
        /*   if ((bin>=num_slices) || (bin<0))  { std::cout<<"bin="<<bin<<"num_slices="<<num_slices<<std::endl;
        throw
        std::runtime_error("something wrong with bining");}*/  

        xkick=dipole_x[bin];
        ykick=dipole_y[bin];
        xkick += -quad_y[bin]*bunch.get_local_particles()[part][0];
        ykick += quad_y[bin]*bunch.get_local_particles()[part][2];     
        zkick = l_monopole[bin];
    
       
        bunch.get_local_particles()[part][1] += wake_factor*xkick;   
        bunch.get_local_particles()[part][3]  += wake_factor*ykick;
        bunch.get_local_particles()[part][5]  += wake_factor*zkick;
      
    }
    
}


void
Impedance::apply(Bunch & bunch, double time_step, Step & step)
{


    double t;
    t = simple_timer_current();
     
    bunch.convert_to_state(Bunch::fixed_t_lab); 
    MArray1d bunchmin, bunchmax;
     
     bunchmin=Diagnostics::calculate_bunchmin(bunch);
     bunchmax=Diagnostics::calculate_bunchmax(bunch);
     double size_z=bunchmax[2]-bunchmin[2];    
     MArray1d  zdensity(boost::extents[z_grid]);
     MArray1d  xmom(boost::extents[z_grid]);
     MArray1d  ymom(boost::extents[z_grid]);
     int lnum_part=bunch.get_local_num();
     MArray1int bin_partition(boost::extents[lnum_part]);
 
     calculate_moments_and_partitions(bunch,zdensity,xmom,ymom,
                      bunchmin[2],size_z, bin_partition);
    
    /*
    std::ofstream file;
    file.open("zdensity.dat");
      for (int i = 0; i < z_grid; ++i){
       file<<i<<"   "<<zdensity[i]<<"   "<<xmom[i]<<"   "<<ymom[i]<<std::endl;
       }  
    file.close();
    abort(); */
     
    
    
     MArray1d  dipole_x(boost::extents[z_grid]);
     MArray1d  dipole_y(boost::extents[z_grid]);
     MArray1d  quad_y(boost::extents[z_grid]);
     MArray1d l_monopole(boost::extents[z_grid]);
    
    
    
   
    double bunchsp=get_bunch_spacing();
    double N_factor=bunch.get_real_num()/bunch.get_total_num();
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta= bunch.get_reference_particle().get_beta();
    double cell_size_z= size_z/double(z_grid);
    
    
  
    std::list<Bunch_means>  stored_bunches = step.get_stored_bunches();    
    get_kicks(z_grid, orbit_length , N_factor, cell_size_z, z_coord, x_wake, y_wake, z_wake, zdensity, xmom, ymom, bunch_spacing,stored_bunches, dipole_x,  dipole_y, quad_y, l_monopole);
   
 


    double w_f=get_wake_factor()*time_step/gamma;  // 1/gamma/beta  is the difference with the old version
                                                        // due to the different units  
 
    impedance_kick(bunch,  w_f, bin_partition, dipole_x, dipole_y, quad_y,  l_monopole);
    t = simple_timer_show(t, "impedance apply");
}

int Impedance::get_z_grid() const { return z_grid;} 
double Impedance::get_orbit_length() const{ return orbit_length;} 
double Impedance::get_wake_factor() const { return wake_factor;} 
double Impedance::get_bunch_spacing() const { return bunch_spacing;}
std::string Impedance::get_pipe_symmetry() const { return pipe_symmetry;}
std::string Impedance::get_wake_file_name() const { return wake_file;}
std::vector<double> Impedance::get_z_coord() const { return z_coord;}
std::vector<double> Impedance::get_x_wake() const { return x_wake;}
std::vector<double> Impedance::get_y_wake() const { return y_wake;}
std::vector<double> Impedance::get_z_wake() const { return z_wake;}
int Impedance::get_nstored_turns() const { return nstored_turns;}

Impedance::~Impedance(){}
