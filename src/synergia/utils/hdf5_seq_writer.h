#ifndef SYNERGIA_UTILS_HDF5_SEQ_WRITER_H
#define SYNERGIA_UTILS_HDF5_SEQ_WRITER_H

#include <vector>
#include <string>

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/multi_array_typedefs.h"

class Hdf5_seq_writer
{
private:

    Hdf5_handler const& file;
    std::string name;
    Commxx const& comm;
    int root;

    int mpi_size;
    int mpi_rank;

    std::vector<hsize_t> fdims;
    std::vector<hsize_t> offset;
    Hdf5_handler dataset;

    bool setup;

public:

    Hdf5_seq_writer(
            Hdf5_handler const& file,
            std::string  const& name,
            Commxx const& comm,
            int root_rank) 
        : file(file)
        , name(name)
        , comm(comm)
        , root(root_rank)
        , mpi_size(comm.size())
        , mpi_rank(comm.rank())
        , fdims()
        , offset()
        , dataset()
        , setup(false)
    { }

    template<typename T>
    void append(T const& data, bool collective)
    {
        auto di = syn::extract_data_info(data);

        if (collective && di.dims.size()==0)
        {
            di.dims = {1};
            //throw std::runtime_error("collective append of scalars not allowed");
        }

        // collect data dims
        auto all_dims0 = syn::collect_dims(
                di.dims, collective, comm, root);

        // offsets for each rank (offsets of dim0 in the combined array)
        std::vector<hsize_t> offsets(mpi_size, 0);
        for (int r=0; r<mpi_size-1; ++r) offsets[r+1] = offsets[r] + all_dims0[r];

        // dim0 of the combined array
        hsize_t dim0 = offsets[mpi_size-1] + all_dims0[mpi_size-1];

        // promote the data dim and set the first dim to 1
        di.dims.resize(di.dims.size() + 1);
        for(int i=di.dims.size()-1; i>0; --i) di.dims[i] = di.dims[i-1];
        di.dims[0] = 1;

        // setup dataset
        if(!setup) 
        {
            setup_dataset(di, offsets, dim0);
            setup = true;
        }
        else
        {
            if (!verify_dims(di, dim0)) throw std::runtime_error(
                    "Hdf5_seq_writer: inconsistent dims in consective write");
        }

        do_append(di, offsets, all_dims0);
    }

    bool verify_dims(syn::data_info_t const& di, hsize_t dim0)
    {
        if (fdims.size()>1 && dim0!=fdims[1]) return false;

        for (int i=2; i<di.dims.size(); ++i)
            if (di.dims[i] != fdims[i]) return false;

        return true;
    }

    void setup_dataset(
            syn::data_info_t& di,
            std::vector<hsize_t> const& offsets,
            hsize_t dim0 )
    {
        const size_t good_chunk_size = 8192; // pulled out of air

        // extended data rank
        auto d_rank = di.dims.size();

        // dims of all data in the current file
        // (num_slabs x size_of_a_slab), or,
        // (num_slabs x all_dim0 x dim1 x dim2 x ...)
        fdims = di.dims;
        if (fdims.size() > 1) fdims[1] = dim0;
        fdims[0] = 0;

        auto max_fdims = fdims;
        auto chunk_fdims = fdims;

        max_fdims[0]   = H5S_UNLIMITED;
        chunk_fdims[0] = (di.data_size < good_chunk_size)
                         ? good_chunk_size/di.data_size : 1;

        // offset of the slab (nslabs, off0, 0, 0, ... )
        offset = std::vector<hsize_t>(d_rank, 0);

        try 
        {
            // try to open the dataset
            dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
            Hdf5_handler fspace = H5Dget_space(dataset);

            // check the data_rank
            if ( H5Sget_simple_extent_ndims(fspace) != d_rank ) 
            {
                throw std::runtime_error(
                        "Hdf5_seq_writer::resumed data has wrong rank");
            }

            std::vector<hsize_t> file_dims(d_rank);
            herr_t res = H5Sget_simple_extent_dims(fspace, &file_dims[0], NULL);
            if (res < 0) throw Hdf5_exception("Error getting dims of the dataspace");

            // check the dims of each rank
            for (int i=1; i<d_rank; ++i)
            {
                if (fdims[i] != file_dims[i]) throw std::runtime_error(
                        "Hdf5_seq_writer::inconsistent data dimensions" );
            }

            fdims[0]  = file_dims[0];
            offset[0] = file_dims[0];
        } 
        catch (Hdf5_exception & e) 
        {
            // dataset not exist, create a new one
            if (dim0 == 0) throw std::runtime_error( 
                    "Hdf5_seq_writer: zero data size encountered");

            Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
            herr_t res = H5Pset_chunk(cparms, d_rank, &chunk_fdims[0]);
            if (res < 0) throw Hdf5_exception();

            Hdf5_handler fspace = H5Screate_simple(d_rank, &fdims[0], &max_fdims[0]);
            dataset = H5Dcreate(file, name.c_str(), di.atomic_type, 
                    fspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
        }
    }


    void do_append(
            syn::data_info_t const& di,
            std::vector<hsize_t> const& offsets,
            std::vector<hsize_t> const& all_dims0 )
    {
        auto dimsm = di.dims;
        if (dimsm.size() > 1) dimsm[1] = all_dims0[mpi_rank];

        // create dataspace for current data block (it looks like max_dims can be null)
        Hdf5_handler mspace = H5Screate_simple(dimsm.size(), dimsm.data(), NULL);

        // increment and extend the dataset to the new size (last_dim+1)
        ++fdims[0];
        herr_t res = H5Dset_extent(dataset, fdims.data());
        if (res < 0) throw Hdf5_exception();

        // select the slab to write
        Hdf5_handler fspace = H5Dget_space(dataset);
        if (offset.size() > 1) offset[1] = offsets[mpi_rank];

        res = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, 
                offset.data(), NULL, dimsm.data(), NULL);

        if (res < 0) throw Hdf5_exception();

        // collective
        Hdf5_handler plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // write
        res = H5Dwrite(dataset, di.atomic_type, mspace, fspace, plist_id, di.ptr);
        if (res < 0) throw Hdf5_exception();

        // increment the offset
        ++offset[0];
    }



};

#endif
