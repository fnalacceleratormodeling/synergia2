#include "synergia/bunch/diagnostics_loss.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/bunch/bunch.h"

Diagnostics_loss::Diagnostics_loss(
        std::string const& filename, 
        std::string const& local_dir )
    : Diagnostics( diag_type, 
                   diag_write_serial, 
                   filename, local_dir )
    , bucket_index(-1)
    , repetition(0)
    , s_ref_particle(0.0)
    , sn_ref_particle(0.0)
    , coords("coords")
{
}  

void
Diagnostics_loss::update(Bunch const& bunch, karray2d_row parts)
{
    auto ref = bunch.get_reference_particle();
    repetition      = ref.get_repetition();
    s_ref_particle  = ref.get_s();
    sn_ref_particle = ref.get_s_n();
    bucket_index    = bunch.get_bucket_index();

    coords = parts;
}


void
Diagnostics_loss::do_write(Bunch const& bunch)
{
    auto & write_helper = get_write_helper(bunch);
    auto const& comm = bunch.get_comm();

    int const mpi_size = comm.size(); 
    int const local_count = coords.extent(0);

    std::vector<int> counts(mpi_size);
    std::vector<int> offsets(mpi_size+1, 0);

    int res = MPI_Gather( &local_count, 1, MPI_INT,
                          counts.data(), 1, MPI_INT,
                          write_helper.get_writer_rank(),
                          comm );

    if (res != MPI_SUCCESS) 
        throw std::runtime_error("MPI error diagnostics_loss::write(): MPI_gather"); 

    for (int i=1; i<mpi_size+1; ++i) 
        offsets[i] = offsets[i-1] + counts[i-1];

    karray2d_row all_coords("all_coords", offsets[mpi_size], 7);

    res = MPI_Gatherv( coords.data(),
                       local_count * 7,
                       MPI_DOUBLE,
                       all_coords.data(),
                       counts.data(),
                       offsets.data(),
                       MPI_DOUBLE,
                       write_helper.get_writer_rank(),
                       comm );

    if (write_helper.write_locally())
    {
        auto & file = write_helper.get_hdf5_file();

        //file.write_serial("bucket_index", bucket_index);
        //file.write_serial("repetition", repetition);
        //file.write_serial("s", s_ref_particle);
        //file.write_serial("s_n", sn_ref_particle);
        file.write_serial("coordinates", all_coords);

        write_helper.finish_write();
    }
}


