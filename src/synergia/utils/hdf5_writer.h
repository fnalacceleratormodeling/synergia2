#ifndef HDF5_WRITER_H_
#define HDF5_WRITER_H_
#include <vector>
#include <string>

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/multi_array_typedefs.h"

class Hdf5_writer
{
private:

    Hdf5_writer() { };

    static void 
    write_impl( Hdf5_handler const& file,
                std::string const& name,
                void const* data_ptr,
                std::vector<hsize_t> const& data_dims,
                std::vector<hsize_t> const& all_dims_0,
                Hdf5_handler const& atomic_type,
                Commxx const& comm );

public:

    // if collective = false, write out a single value to the file (value from the root
    // rank)
    //
    // if collective = true, do a manual reduction (only in serial hdf5 ), and write 
    // out the data to hdf5 file. scalar data will be extended a 1-d array, vector data 
    // of 1-d or higher dim will get reduced on the first dimension. the extent of the 
    // first dim can be different, but all higher dims must be exactly the same across 
    // all ranks
    //
    // e.g.: 
    //
    //   r0->write(10, false),      r1->write(20, false)  => file0: 10
    //   r0->write(10, true),       r1->write(20, true)   => file0: [10, 20]
    //   r0->write([10, 11], true), r1->write([20], true) => file0: [10, 11, 20]

    template<class T>
    static void 
    write( Hdf5_handler const& file, 
           std::string const& name, 
           T const& data, 
           bool collective, 
           Commxx const& comm,
           int root_rank );

    template<class T>
    static void 
    write( Hdf5_handler const& file, 
           std::string const& name,
           T const* data, size_t len, 
           bool collective, 
           Commxx const& comm,
           int root_rank );
};

template<class T>
inline void 
Hdf5_writer::write( 
        Hdf5_handler const& file, 
        std::string const& name, 
        T const& data, 
        bool collective, 
        Commxx const& comm,
        int root_rank )
{
    auto di = syn::extract_data_info(data);

    // promote the dimensionality if collective write a scaler
    if (collective && di.dims.size()==0) 
        di.dims = {1};

    auto all_dim0 = syn::collect_dims(
            di.dims, collective, comm, root_rank);

    write_impl(file, name, di.ptr, di.dims, 
            all_dim0, di.atomic_type, comm);
}

template<class T>
inline void 
Hdf5_writer::write( 
        Hdf5_handler const& file, 
        std::string const& name,
        T const* data, size_t len, 
        bool collective, 
        Commxx const& comm,
        int root_rank )
{
    auto data_ptr  = data;
    auto data_dims = std::vector<hsize_t>({len});
    auto atomic    = hdf5_atomic_data_type<T>();

    auto all_dim0 = syn::collect_dims(
            data_dims, collective, comm, root_rank);

    write_impl(file, name, data_ptr, data_dims, 
            all_dim0, atomic, comm);
}

inline void 
Hdf5_writer::write_impl( 
        Hdf5_handler const& file,
        std::string const& name,
        void const* data_ptr,
        std::vector<hsize_t> const& data_dims,
        std::vector<hsize_t> const& all_dims_0,
        Hdf5_handler const& atomic_type,
        Commxx const& comm )
{
    int mpi_size = comm.size();
    int mpi_rank = comm.rank();

    // data rank
    auto data_rank = data_dims.size();

    // offsets for each rank (offsets of dim0 in the combined array)
    std::vector<hsize_t> offsets(mpi_size, 0);
    for (int r=0; r<mpi_size-1; ++r) offsets[r+1] = offsets[r] + all_dims_0[r];

    // dim0 of the combined array
    hsize_t dim0 = offsets[mpi_size-1] + all_dims_0[mpi_size-1];

    // dims for the combined array
    std::vector<hsize_t> dimsf(data_dims);
    if (dimsf.size()) dimsf[0] = dim0;

#ifdef USE_PARALLEL_HDF5

    if (!file.valid())
        throw std::runtime_error("invalid file handler");

    // dataset
    Hdf5_handler filespace = H5Screate_simple(data_rank, dimsf.data(), NULL);
    Hdf5_handler dset = H5Dcreate(file, name.c_str(), atomic_type,
            filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // only create the dataset, but do not initiate the write 
    // if the total size is 0
    if (dimsf.size() && dimsf[0] == 0) return;

    // local dims(counts)
    auto dimsm = data_dims;
    if (dimsm.size()) dimsm[0] = all_dims_0[mpi_rank];

    // dataspace
    Hdf5_handler memspace = H5Screate_simple(data_rank, dimsm.data(), NULL);
    Hdf5_handler filespace2 = H5Dget_space(dset);

    // select hyperslab only for non-scalars
    if (data_rank)
    {
        auto offset = std::vector<hsize_t>(data_rank, 0);
        offset[0] = offsets[mpi_rank];

        herr_t res = H5Sselect_hyperslab(filespace2, H5S_SELECT_SET, 
                offset.data(), NULL, dimsm.data(), NULL);
        
        if (res < 0) throw Hdf5_exception();
    }

    // collective write
    Hdf5_handler plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // write
    herr_t res = H5Dwrite(dset, atomic_type, memspace, filespace2, plist_id, data_ptr);
    if (res < 0) throw Hdf5_exception();

#else

    // TODO:
#if 0
    if (comm.rank() == write_rank)
    {
        if (!file.valid())
            throw std::runtime_error("invalid file handler");

        Hdf5_handler filespace = H5Screate_simple(data_rank, dimsf, NULL);
        Hdf5_handler dset = H5Dcreate(file, name.c_str(), atomic_type,
                filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for(r=0; i<comm.size(); ++r)
        {
            dimsm = all_dims[r];
            offset = all_offsets[r];

            if (dimsm == 0) continue;

            Hdf5_handler memspace = H5Screate_simple(data_rank, dimsm, NULL);
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dimsm, NULL);

            if (r!= comm.rank())
            {
                mpi_recv(recv_data, r);
            }
            else
            {
                recv_data = transpose(data);
            }

            H5Dwrite(dset, atomic, memspace, filespace, default, recv_data);

        }
    }
    else
    {
        if (collective)
        {
            data_buf = transpose(data);
            mpi_send(data_buf, write_rank);
        }
    }
#endif

#endif
}


#endif /* HDF5_WRITER_H_ */
