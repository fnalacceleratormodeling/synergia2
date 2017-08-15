#include "diagnostics_deltae_particles.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <boost/filesystem.hpp>

const char Diagnostics_deltae_particles::name[] = "diagnostics_deltae_particles";

Diagnostics_deltae_particles::Diagnostics_deltae_particles(std::string const& filename,
                                     int begin_turn, int end_turn,  int average_bin,
                                     double beta_x, double alpha_x, double beta_y, double alpha_y,                    
                                     double beta_z, double alpha_z, std::string const& local_dir):
Diagnostics(Diagnostics_deltae_particles::name, filename, local_dir), begin_turn(begin_turn), end_turn(end_turn),
average_bin(average_bin), beta_x(beta_x), alpha_x(alpha_x), beta_y(beta_y), alpha_y(alpha_y), 
beta_z(beta_z), alpha_z(alpha_z), have_writers(false), have_initial_deltae(false)
{
    if ((end_turn-begin_turn)<average_bin){
      throw std::runtime_error("diagnostics_deltae_particles: (end_turn-begin_turn) is smaller than average_bin)");
    }
    
}


Diagnostics_deltae_particles::Diagnostics_deltae_particles()
{
}

bool
Diagnostics_deltae_particles::is_serial() const
{
    return false;
}

void
Diagnostics_deltae_particles::update()
{

     int local_num=0; 
     if ((this->have_bunch()) && (!have_initial_deltae )){
        if (get_bunch().get_comm().has_this_rank()){
          local_num = get_bunch().get_local_num(); 
        }
        
        local_deltae = new MArray2d(boost::extents[local_num][4]);
        for (int part = 0; part < local_num; ++part){
              (*local_deltae)[part][0]=0.;
              (*local_deltae)[part][1]=0.;
              (*local_deltae)[part][2]=0.;
              (*local_deltae)[part][3]=0.;
           }
      turn=0; 
      have_initial_deltae=true; 
     }
     
   
    if (get_bunch().get_comm().has_this_rank()){
        turn +=1;
        get_bunch().convert_to_state(get_bunch().fixed_z_lab);
        int local_num = get_bunch().get_local_num();
        Const_MArray2d_ref particles(get_bunch().get_local_particles());
         if ((turn>=begin_turn) && (turn<end_turn)){
      
            for (int part = 0; part < local_num; ++part) {
                (*local_deltae)[part][1] +=(particles[part][Bunch::x]*particles[part][Bunch::x]
                        +(alpha_x*particles[part][Bunch::x]+beta_x*particles[part][Bunch::xp])
                          *(alpha_x*particles[part][Bunch::x]+beta_x*particles[part][Bunch::xp]))/(2.*beta_x*(end_turn-begin_turn));   
                 (*local_deltae)[part][2] +=  (particles[part][Bunch::y]*particles[part][Bunch::y]
                          +(alpha_y*particles[part][Bunch::y]+beta_y*particles[part][Bunch::yp])
                            *(alpha_y*particles[part][Bunch::y]+beta_y*particles[part][Bunch::yp]))/(2.*beta_y*(end_turn-begin_turn));    
              
                (*local_deltae)[part][3] +=  (particles[part][Bunch::z]*particles[part][Bunch::z]
                          +(alpha_z*particles[part][Bunch::z]+beta_z*particles[part][Bunch::zp])
                            *(alpha_z*particles[part][Bunch::z]+beta_z*particles[part][Bunch::zp]))/(2.*fabs(beta_z)*(end_turn-begin_turn));            
                            
            }
         }
        
        if ((turn>=begin_turn) && (turn<begin_turn+average_bin)) {
            for (int part = 0; part < local_num; ++part) {
              (*local_deltae)[part][0] -=
              (particles[part][Bunch::x]*particles[part][Bunch::x]
                +(alpha_x*particles[part][Bunch::x]+beta_x*particles[part][Bunch::xp])
                     *(alpha_x*particles[part][Bunch::x]+beta_x*particles[part][Bunch::xp]))/(2.*beta_x*beta_x*average_bin)
              +(particles[part][Bunch::y]*particles[part][Bunch::y]
                  +(alpha_y*particles[part][Bunch::y]+beta_y*particles[part][Bunch::yp])
                      *(alpha_y*particles[part][Bunch::y]+beta_y*particles[part][Bunch::yp]))/(2.*beta_y*beta_y*average_bin);  
                           
           }      
        }
        else if ((turn>=end_turn-average_bin) && (turn<end_turn)){
            for (int part = 0; part < local_num; ++part) {
              (*local_deltae)[part][0] += 
                (particles[part][Bunch::x]*particles[part][Bunch::x]
                +(alpha_x*particles[part][Bunch::x]+beta_x*particles[part][Bunch::xp])
                     *(alpha_x*particles[part][Bunch::x]+beta_x*particles[part][Bunch::xp]))/(2.*beta_x*beta_x*average_bin)
              +(particles[part][Bunch::y]*particles[part][Bunch::y]
                  +(alpha_y*particles[part][Bunch::y]+beta_y*particles[part][Bunch::yp])
                      *(alpha_y*particles[part][Bunch::y]+beta_y*particles[part][Bunch::yp]))/(2.*beta_y*beta_y*average_bin); 
             
           }  
        }
    }
}

void
Diagnostics_deltae_particles::send_local_particles()
{
    if (get_bunch().get_comm().has_this_rank()){
        int local_num = get_bunch().get_local_num();
        void * send_buffer = (void*) (*local_deltae).origin();          
        int status;
        int message_size =4*local_num;
        int receiver = get_write_helper().get_writer_rank();
        int rank = get_bunch().get_comm().get_rank();
        MPI_Comm comm = get_bunch().get_comm().get();
        status = MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
            comm);
        if (status != MPI_SUCCESS) {
            throw std::runtime_error(
                "Diagnostics_deltae_particles::send_local_particles: MPI_Send failed.");
        }
    }
}


void
write_selected_particles(Hdf5_chunked_array2d_writer & writer,
        MArray2d_ref const & particles, int local_num)
{   
        writer.write_chunk(
                particles[boost::indices[range(0, local_num)][range()]]);  
}


void
Diagnostics_deltae_particles::receive_other_local_particles(
        std::vector<int > const& local_nums, Hdf5_file_sptr file_sptr)
{
    if (get_bunch().get_comm().has_this_rank()){
     int myrank = get_bunch().get_comm().get_rank();
     int size = get_bunch().get_comm().get_size();
     Hdf5_chunked_array2d_writer
         writer_particles(
             &(file_sptr->get_h5file()),
             "deltae_and_Jave",*local_deltae);
            
     for (int rank = 0; rank < size; ++rank) {
         int local_num = local_nums[rank];
         if (rank == myrank) {
         write_selected_particles(writer_particles,
             *local_deltae, local_num );
         } else {
         MPI_Status status;
         Raw_MArray2d received(boost::extents[local_num][4]);
         int message_size = 4*local_num;
         MPI_Comm comm = get_bunch().get_comm().get();
         int error = MPI_Recv((void*) received.m.origin(), message_size,
             MPI_DOUBLE, rank, rank, comm, &status);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                "Diagnostics_deltae_particles::receive_other_local_particles: MPI_Recv failed.");
        }
         write_selected_particles(writer_particles, received.m, local_num);
         }
     }
    }
}

void
Diagnostics_deltae_particles::write()
{
  if (turn==end_turn+1){
      if (get_bunch().get_comm().has_this_rank()){
        get_bunch().convert_to_state(get_bunch().fixed_z_lab);
        
        MPI_Comm comm = get_bunch().get_comm().get();
        int local_num = get_bunch().get_local_num();
        int num_procs = get_bunch().get_comm().get_size();
        std::vector<int > local_nums(num_procs);
        void * local_nums_buf = (void *) &local_nums[0];
        int root = get_write_helper().get_writer_rank();
        int status;
        status = MPI_Gather((void*) &local_num, 1, MPI_INT, local_nums_buf, 1,
            MPI_INT, root, comm);
        if (status != MPI_SUCCESS) {
            throw std::runtime_error(
                "Diagnostics_particles::write: MPI_Gather failed.");
        }
        

      if (get_write_helper().write_locally()) {
            Hdf5_file_sptr file_sptr = get_write_helper().get_hdf5_file_sptr();
            receive_other_local_particles(local_nums, file_sptr);
            get_write_helper().finish_write();
        } else {
            send_local_particles();
            }
        
          
      }
  }
  
}


std::string
Diagnostics_deltae_particles::get_deltae_serialization_path() const
{
    
    std::stringstream sstream;
    sstream << "deltae_";
    sstream <<get_bunch().get_bucket_index();
    sstream << ".h5";
    return get_serialization_path(sstream.str());
}


    
template<class Archive>
    void
    Diagnostics_deltae_particles::save(Archive & ar, const unsigned int version) const
    {     
        ar << BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
            << BOOST_SERIALIZATION_NVP(have_writers)
            << BOOST_SERIALIZATION_NVP(turn)
            << BOOST_SERIALIZATION_NVP(begin_turn)
            << BOOST_SERIALIZATION_NVP(end_turn)
            << BOOST_SERIALIZATION_NVP(average_bin)
            << BOOST_SERIALIZATION_NVP(beta_x)
            << BOOST_SERIALIZATION_NVP(alpha_x)
            << BOOST_SERIALIZATION_NVP(beta_y)
            << BOOST_SERIALIZATION_NVP(alpha_y)
            << BOOST_SERIALIZATION_NVP(beta_z)
            << BOOST_SERIALIZATION_NVP(alpha_z)
            << BOOST_SERIALIZATION_NVP(have_initial_deltae);                      
        if (get_bunch().get_comm().has_this_rank()){
            int attempts=0;
            bool fail=true;           
            while ((attempts<5) && fail){

                try {
                boost::filesystem::remove(get_deltae_serialization_path());
                    Hdf5_file file(get_deltae_serialization_path(), Hdf5_file::truncate);                   
                     file.write(*local_deltae, "local_deltae");
                     file.close();
                     fail=false;
                }
                catch(H5::Exception& he) {
                    ++attempts;
                    fail=true;
                    std::cout<<"diagnostics_deltae: H5 Exception thrown, attempts number="
                        <<attempts<<" on rank="<<Commxx().get_rank()<<std::endl;
                    sleep(3);
                }
            }                   
        }
            
    }
    
  template<class Archive>
    void
    Diagnostics_deltae_particles::load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
            >> BOOST_SERIALIZATION_NVP(have_writers)
            >> BOOST_SERIALIZATION_NVP(turn)
            >> BOOST_SERIALIZATION_NVP(begin_turn)
            >> BOOST_SERIALIZATION_NVP(end_turn)
            >> BOOST_SERIALIZATION_NVP(average_bin)
            >> BOOST_SERIALIZATION_NVP(beta_x)
            >> BOOST_SERIALIZATION_NVP(alpha_x)
            >> BOOST_SERIALIZATION_NVP(beta_y)
            >> BOOST_SERIALIZATION_NVP(alpha_y)
            >> BOOST_SERIALIZATION_NVP(beta_z)
            >> BOOST_SERIALIZATION_NVP(alpha_z)
            >> BOOST_SERIALIZATION_NVP(have_initial_deltae);                     
       if (get_bunch().get_comm().has_this_rank()) {
            Hdf5_file file(get_deltae_serialization_path(),
                    Hdf5_file::read_only);             
            local_deltae
                    = new MArray2d(file.read<MArray2d > ("local_deltae"));
        } else {
           int local_num = 0;
           local_deltae = new MArray2d(boost::extents[local_num][4]);
        }   
        
    }
    

template
void
Diagnostics_deltae_particles::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_deltae_particles::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_deltae_particles::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_deltae_particles::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_deltae_particles::~Diagnostics_deltae_particles()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_deltae_particles)
