#include "diagnostics_phase_space_density.h"
#include "synergia/utils/simple_timer.h"

const char Diagnostics_phase_space_density::name[] = "diagnostics_phase_space_density";


Diagnostics_phase_space_density::Diagnostics_phase_space_density(std::string const& filename,
                  int grid_z, int grid_zp,double  z_nsigma, double zp_nsigma,  std::string const& local_dir):
  Diagnostics_phase_space_density:: Diagnostics(Diagnostics_phase_space_density::name, filename,
                local_dir),  grid_z(grid_z), grid_zp(grid_zp), z_nsigma(z_nsigma), 
                zp_nsigma(zp_nsigma), have_grid(false), have_writers(false), 
                writer_s_n(0), writer_repetition(0), writer_s(0), density(boost::extents[grid_z][grid_zp]), 
                writer_density(0), xdensity(boost::extents[grid_z][grid_zp]), writer_xdensity(0)
               // ydensity(boost::extents[grid_z][grid_zp]), writer_ydensity(0)
{
}

Diagnostics_phase_space_density::Diagnostics_phase_space_density(): have_writers(false)
{
}


void
Diagnostics_phase_space_density::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
         file_sptr->write(grid_z, "grid_z");
         file_sptr->write(grid_zp, "grid_zp");
         file_sptr->write(z_nsigma, "z_nsigma");
         file_sptr->write(zp_nsigma, "zp_nsigma");
         if (have_grid) {
           file_sptr->write(grid_zrms, "grid_zrms");
           file_sptr->write(grid_zprms, "grid_zprms");
         }
         else{
           throw std::runtime_error(
                "update shoule be called before init_writers");
         }
         

        writer_s_n = new Hdf5_serial_writer<double >(file_sptr, "s_n");
        writer_repetition = new Hdf5_serial_writer<int >(file_sptr,
                "repetition");
        writer_s = new Hdf5_serial_writer<double >(file_sptr,
                "s");
        writer_density = new Hdf5_serial_writer<MArray2d_ref  >(file_sptr,
                "density");  
        writer_xdensity = new Hdf5_serial_writer<MArray2d_ref  >(file_sptr,
                "xdensity");   
       // writer_ydensity = new Hdf5_serial_writer<MArray2d_ref  >(file_sptr,
       //         "ydensity");        
        have_writers = true;
    }
}




bool
Diagnostics_phase_space_density::is_serial() const
{
    return true;
}



void
Diagnostics_phase_space_density::update()
{ 
  if (get_bunch().get_comm().has_this_rank()){
    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
    s_n = get_bunch().get_reference_particle().get_s_n();
    repetition = get_bunch().get_reference_particle().get_repetition();
    s= get_bunch().get_reference_particle().get_s();
    if (!have_grid){
      MArray1d mean = Core_diagnostics::calculate_mean(get_bunch());
      MArray1d std = Core_diagnostics::calculate_std(get_bunch(), mean);
     // std::cout<<" stdz="<<std[4]<<" stdzp="<<std[5]<<std::endl;
      grid_zrms=std[4];
      grid_zprms=std[5];
      z_range=z_nsigma*grid_zrms;
      z_cell= z_range/grid_z;
      zp_range=zp_nsigma*grid_zprms;
      zp_cell= zp_range/grid_zp;
    //  std::cout<<" grid_zrms="<<grid_zrms<<" grid_zprms="<<grid_zprms<<std::endl;
      
      have_grid=true;
    }
    deposit_densities();        
  }
}

void
Diagnostics_phase_space_density::write()
{
   if (get_bunch().get_comm().has_this_rank()){
     if (get_write_helper().write_locally()) {
        init_writers(get_write_helper().get_hdf5_file_sptr());      
        writer_s_n->append(s_n);
        writer_repetition->append(repetition);
        writer_s->append(s);
        writer_density->append(density);
        writer_xdensity->append(xdensity);
       // writer_ydensity->append(ydensity);
        get_write_helper().finish_write();
    }     
   }
}  

void
Diagnostics_phase_space_density::deposit_densities()
{ 
  
    for (int i = 0; i < grid_z; ++i) {
      for (int j = 0; j < grid_zp; ++j) {
          density[i][j]=0.;
          xdensity[i][j]=0.;
         // ydensity[i][j]=0.;
      }        
    }
    
    double normalization=1./get_bunch().get_total_num();
    MArray1d mean = Core_diagnostics::calculate_mean(get_bunch());
    double x_mean=mean[0];
  //  double y_mean=mean[2];
    double z_mean=mean[4];
    double zp_mean=mean[5];
    Const_MArray2d_ref particles(get_bunch().get_local_particles());
    for (int part = 0; part < get_bunch().get_local_num(); ++part) {
        double z_dist=particles[part][4]-z_mean+0.5*z_range;
        double  zp_dist=particles[part][5]-zp_mean+0.5*zp_range;
         if ((z_dist>0) && (z_dist<z_range) 
                 && (zp_dist>0) && (zp_dist<zp_range)){
             
             int iz= int(floor(z_dist/z_cell-0.5));
             int izp= int(floor(zp_dist/zp_cell-0.5));
             double off_iz=z_dist/z_cell-iz-0.5;
             double off_izp=zp_dist/zp_cell-izp-0.5;
             for (int ih = 0; ih < 2; ++ih) {
                int cellx = iz + ih;
                for (int iv = 0; iv < 2; ++iv) { 
                   int celly = izp + iv;
                    if ((cellx>=0) && (celly>=0) && (cellx<grid_z) && (celly<grid_zp)){
                      double weight=(1 - ih - (1 - 2 * ih)* off_iz) * (1 - iv - (1 - 2 * iv) * off_izp);
                      density[cellx][celly] += weight*normalization;
                      xdensity[cellx][celly] += weight*normalization*(particles[part][0]-x_mean);
                     // ydensity[cellx][celly] += weight*normalization*(particles[part][2]-y_mean);
                    }
                  
                }
             }          
         }       
    }// for part
    
   
     int error = MPI_Allreduce(MPI_IN_PLACE, (void*)  density.origin(),
                               density.num_elements(), MPI_DOUBLE, MPI_SUM, get_bunch().get_comm().get());
                               
     if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in diagnostics_phase_space_density: MPI_Allreduce in deposit_density1");
    }
    
     error = MPI_Allreduce(MPI_IN_PLACE, (void*)  xdensity.origin(),
                               xdensity.num_elements(), MPI_DOUBLE, MPI_SUM, get_bunch().get_comm().get());
    
    
     if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in diagnostics_phase_space_density: MPI_Allreduce in deposit_xdensity");
    }
    
   //  error = MPI_Allreduce(MPI_IN_PLACE, (void*)  ydensity.origin(),
   //                            ydensity.num_elements(), MPI_DOUBLE, MPI_SUM, get_bunch().get_comm().get());
    
    
   //  if (error != MPI_SUCCESS) {
   //     throw std::runtime_error(
   //             "MPI error in diagnostics_phase_space_density: MPI_Allreduce in deposit_ydensity");
  //  }
    
    
}

template<class Archive>
    void
    Diagnostics_phase_space_density::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
        ar & BOOST_SERIALIZATION_NVP(grid_z);                       
        ar & BOOST_SERIALIZATION_NVP(grid_zp);                      
        ar & BOOST_SERIALIZATION_NVP(z_nsigma);                     
        ar & BOOST_SERIALIZATION_NVP(zp_nsigma);                    
        ar & BOOST_SERIALIZATION_NVP(grid_zrms);                    
        ar & BOOST_SERIALIZATION_NVP(grid_zprms);                   
        ar & BOOST_SERIALIZATION_NVP(z_range);                      
        ar & BOOST_SERIALIZATION_NVP( zp_range);                    
        ar & BOOST_SERIALIZATION_NVP(z_cell);                       
        ar & BOOST_SERIALIZATION_NVP(zp_cell);                      
        ar & BOOST_SERIALIZATION_NVP(have_grid);                    
        ar & BOOST_SERIALIZATION_NVP(have_writers);                 
        ar & BOOST_SERIALIZATION_NVP(s_n);                          
        ar & BOOST_SERIALIZATION_NVP(writer_s_n);                   
        ar & BOOST_SERIALIZATION_NVP(repetition);                   
        ar & BOOST_SERIALIZATION_NVP(writer_repetition);            
        ar & BOOST_SERIALIZATION_NVP(s);                            
        ar & BOOST_SERIALIZATION_NVP(writer_s);                     
        ar & BOOST_SERIALIZATION_NVP(density);                      
        ar & BOOST_SERIALIZATION_NVP(writer_density);               
        ar & BOOST_SERIALIZATION_NVP(xdensity);                     
        ar & BOOST_SERIALIZATION_NVP(writer_xdensity); 
      //  ar & BOOST_SERIALIZATION_NVP(ydensity);                     
       // ar & BOOST_SERIALIZATION_NVP(writer_ydensity);
                                                                    
    }
                                                                  
template                                                          
void                                                              
Diagnostics_phase_space_density::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);
                                                                  
template                                                          
void                                                              
Diagnostics_phase_space_density::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);
                                                                  
template                                                          
void                                                              
Diagnostics_phase_space_density::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
                                                                  
template                                                          
void                                                              
Diagnostics_phase_space_density::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);



Diagnostics_phase_space_density::~Diagnostics_phase_space_density()
{
  if (have_writers) {
        delete writer_s;
        delete writer_repetition;
        delete writer_s_n;
        delete writer_density;
        delete writer_xdensity;
       // delete writer_ydensity;
       
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_phase_space_density)