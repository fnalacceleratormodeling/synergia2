#include "diagnostics_loss.h"
#include "synergia/foundation/math_constants.h"

Diagnostics_loss::Diagnostics_loss(
        std::string const& filename, 
        std::string const& local_dir )
    : Diagnostics( diag_type, 
                   diag_write_serial, 
                   filename, local_dir )
{
}  


#if 0
void
Diagnostics_loss::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
           writer_bucket_index = new Hdf5_serial_writer<int  > (file_sptr,  "bucket_index");
           writer_repetition = new Hdf5_serial_writer<int  > (file_sptr,  "repetition");
           writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");    
           writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
           writer_coords = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "coordinates");
           have_writers = true;
    }
}
#endif





void
Diagnostics_loss::update(karray2d_row parts)
{
    Bunch const& b = get_bunch();

#if 0
    bucket_index.push_back(b.get_bucket_index());
    repetition.push_back(b.get_rep());
    s_ref_particle.push_back(b.get_s());
    sn_ref_particle.push_back(sn);
#endif
    coords = parts;
}


void
Diagnostics_loss::do_write()
{
#if 0
    if (!have_bunch()) { throw std::runtime_error(
          "MPI error diagnostics_apertures_loss:write  the diagnostics does not have a bunch_sptr"); 
    }else{

        Commxx_sptr comm_sptr=get_bunch().get_comm_sptr();
        int const mpi_size=comm_sptr->get_size(); 
        int const rank=comm_sptr->get_rank();
        int local_size=bucket_index.size();
        std::vector<int> counts(mpi_size,1);
        std::vector<int> offsets(mpi_size,0);

        for (int r = 1; r < mpi_size; ++r){
              offsets[r] =offsets[r-1]+counts[r-1] ;
        }
        
        MArray1i data_sizes(boost::extents[mpi_size]);       
        int error = MPI_Allgatherv((void *) &local_size, counts[rank], MPI_INT,
            reinterpret_cast<void*>(data_sizes.origin()), &counts[0], &offsets[0], MPI_INT, comm_sptr->get());
        if (error != MPI_SUCCESS) {
                throw std::runtime_error(
              "MPI error diagnostics_apertures_write:MPI_Allgatherv 1"); 
        }       
          
        offsets[0]=0; 
        for (int r = 1; r < mpi_size; ++r){
              offsets[r] =offsets[r-1]+data_sizes[r-1] ;
        }


        int nlost= offsets[mpi_size-1]+data_sizes[mpi_size-1];
        MArray1i bucket_index_array(boost::extents[nlost]);
        error = MPI_Allgatherv((void *) &bucket_index[0], local_size, MPI_INT,
            reinterpret_cast<void*>(bucket_index_array.origin()), &data_sizes[0], &offsets[0], MPI_INT, comm_sptr->get());
        if (error != MPI_SUCCESS) {
                throw std::runtime_error(
              "MPI error diagnostics_apertures_write:MPI_Allgatherv bucket_index"); 
        }       
          
        MArray1i repetition_array(boost::extents[nlost]);
        error = MPI_Allgatherv((void *) &repetition[0], local_size, MPI_INT,
            reinterpret_cast<void*>(repetition_array.origin()), &data_sizes[0], &offsets[0], MPI_INT, comm_sptr->get());
        if (error != MPI_SUCCESS) {
                throw std::runtime_error(
              "MPI error diagnostics_apertures_write:MPI_Allgatherv repetition"); 
        }  
        
        MArray1d s_ref_particle_array(boost::extents[nlost]);
        error = MPI_Allgatherv((void *) &s_ref_particle[0], local_size, MPI_DOUBLE,
            reinterpret_cast<void*>(s_ref_particle_array.origin()), &data_sizes[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
        if (error != MPI_SUCCESS) {
                throw std::runtime_error(
              "MPI error diagnostics_apertures_write:MPI_Allgatherv s_ref_particle"); 
        }  
        
        MArray1d sn_ref_particle_array(boost::extents[nlost]);
        error = MPI_Allgatherv((void *) &sn_ref_particle[0], local_size, MPI_DOUBLE,
            reinterpret_cast<void*>(sn_ref_particle_array.origin()), &data_sizes[0], &offsets[0], MPI_DOUBLE, comm_sptr->get());
        if (error != MPI_SUCCESS) {
                throw std::runtime_error(
              "MPI error diagnostics_apertures_write:MPI_Allgatherv sn_ref_particle"); 
        }  
        
        MArray2d local_coords(boost::extents[local_size][7]);
        for (int i = 0; i < local_size; ++i){
          for (int j = 0; j < 7; ++j){
              local_coords[i][j]=coords[i][j];
          }
        }
        std::vector<int> coord_counts(mpi_size);
        std::vector<int> coord_offsets(mpi_size);
        for (int r = 0; r < mpi_size; ++r){
              coord_offsets[r] =offsets[r]*7;
              coord_counts[r]=data_sizes[r]*7;
        }
        
        
        MArray2d coords_array(boost::extents[nlost][7]);
        error = MPI_Allgatherv((void *) local_coords.data(), 7*local_size, MPI_DOUBLE,
            reinterpret_cast<void*>(coords_array.origin()),  & coord_counts[0], &coord_offsets[0], MPI_DOUBLE, comm_sptr->get());
        if (error != MPI_SUCCESS) {
                throw std::runtime_error(
              "MPI error diagnostics_apertures_write:MPI_Allgatherv coords"); 
        }  

      
       
            if (get_bunch().get_comm().has_this_rank()){  
                if  (nlost>0) {
                    if (get_write_helper().write_locally()) {                  
                            init_writers(get_write_helper().get_hdf5_file_sptr());
                            for (int part = 0; part < nlost; ++part){
                                writer_bucket_index->append(bucket_index_array[part]);
                                writer_repetition->append(repetition_array[part]);  
                                writer_s->append(s_ref_particle_array[part]); 
                                writer_s_n->append(sn_ref_particle_array[part]);
                                MArray1d coords_view(coords_array[part][boost::indices[range(0, 7)]]);              
                                writer_coords->append(coords_view);
                            }
                            get_write_helper().finish_write();
                    }
               }
           }
                          
    
        bucket_index.clear();
        repetition.clear();
        s_ref_particle.clear();
        sn_ref_particle.clear();
        coords.clear();
    }            
#endif
}


