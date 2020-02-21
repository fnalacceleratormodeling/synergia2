#ifndef HDF5_READER_H
#define HDF5_READER_H

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/multi_array_typedefs.h"

class Hdf5_reader
{
public:

    Hdf5_reader() = delete;

    // if T is a scalar type, raed in a single value
    // if T is kokkos::view, read in the array data
    //
    // notes:
    // 1. the array container must be in c ordering
    //
    // 2. the read method does a dimension check. e.g.,
    //    a 3d dataset can be read back only via a
    //    3d container (kokkos::view) as the template
    //    parameter
    //
    // 3. the returned container will have a shape
    //    suitable for holding the entire dataset
    //
    // 4. no collective read allowed for now. the same 
    //    data object containing the entire array will
    //    be returned on every rank
    //
    template<class T>
    static T
    read( Hdf5_handler const& file,
          std::string const& name,
          Commxx const& comm,
          int root_rank )
    {
        T retval;
        auto di = syn::extract_data_info(retval);

        // collect dim0
        auto all_dim0 = syn::collect_dims(
                di.dims, false, comm, root_rank );

        read_verify_and_set_data_dims(
                file, name, all_dim0, di,
                false, comm, root_rank );

        syn::resize_data_obj(retval, di);

        read_impl(file, name, all_dim0, di, comm);

        return retval;
    }

    template<class T>
    static T
    read( Hdf5_handler const& file,
          std::string const& name,
          size_t len,
          Commxx const& comm,
          int root_rank )
    {
        T retval;
        auto di = syn::extract_data_info(retval);

        // collective read from an array into a
        // scalar doesnt make much sense
        if (di.dims.size()==0) throw std::runtime_error(
                "hdf5_reader::read() read from an array into a scalar");

        // set the first dim, rests are still 0s
        di.dims[0] = len;

        // collect
        auto all_dim0 = syn::collect_dims(
                di.dims, true, comm, root_rank );

        read_verify_and_set_data_dims(
                file, name, all_dim0, di,
                true, comm, root_rank );

        syn::resize_data_obj(retval, di);

        read_impl(file, name, all_dim0, di, comm);

        return retval;
    }

    // read 1d array
    template<class T>
    static void
    read( Hdf5_handler const& file,
          std::string const& name,
          T * data, 
          Commxx const& comm,
          int root_rank )
    {
        throw std::runtime_error("hdf5_reader::read_3() not implemented");
    }
 
    // collective read a 1d array
    template<class T>
    static void
    read( Hdf5_handler const& file,
          std::string const& name,
          T * data, 
          size_t len,
          Commxx const& comm,
          int root_rank )
    {
        syn::data_info_t di{ 
            data, {len}, hdf5_atomic_data_type<T>() };

        auto all_dim0 = syn::collect_dims(
                di.dims, true, comm, root_rank );

        read_verify_and_set_data_dims(
                file, name, all_dim0, di,
                true, comm, root_rank );

        read_impl(file, name, all_dim0, di, comm);
    }

    static void 
    read_verify_and_set_data_dims(
            Hdf5_handler const& file,
            std::string  const& name,
            std::vector<hsize_t>& all_dims_0,
            syn::data_info_t& di,
            bool collective,
            Commxx const& comm,
            int root_rank )
    {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        // data rank
        auto data_rank = di.dims.size();

        // offsets for each rank (offsets of dim0 in the combined array)
        std::vector<hsize_t> offsets(mpi_size, 0);
        for (int r=0; r<mpi_size-1; ++r) offsets[r+1] = offsets[r] + all_dims_0[r];

        // dim0 of the combined array
        hsize_t dim0 = offsets[mpi_size-1] + all_dims_0[mpi_size-1];

        // get and verify done only on the root_rank
        if (mpi_rank == root_rank)
        {
            if (!file.valid())
                throw std::runtime_error("invalid file handler");

            Hdf5_handler dset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
            Hdf5_handler filespace = H5Dget_space(dset);

            // get and check the datatype
            Hdf5_handler type = H5Dget_type(dset);
            if (H5Tequal(type, di.atomic_type) <= 0)
                throw std::runtime_error("wrong atomic data type");

            // get and check the data rank
            int rank = H5Sget_simple_extent_ndims(filespace);
            if (rank < 0) throw Hdf5_exception("error at getting rank");

            if (rank != data_rank) 
                throw std::runtime_error("Hdf5_reader: wrong rank");

            // get, check, and broadcast the dims for non-scalar
            if (rank)
            {
                std::vector<hsize_t> dims(rank);
                herr_t res = H5Sget_simple_extent_dims(filespace, dims.data(), NULL);
                if (res < 0) throw Hdf5_exception("error at getting dims");

                // only the collective read of non-scalar needs the dim check
                if (collective && (dims[0] != dim0))
                {
                    throw std::runtime_error(
                            "Hdf5_reader: combined dim[0] is inconsistent with "
                            "the dataset from file" );

                }

                di.dims = dims;
            }
        }

        // broadcast the dims to all ranks
        if (data_rank)
            MPI_Bcast(di.dims.data(), data_rank, MPI_UINT64_T, root_rank, comm);

        // in collective vector read, all_dim0 is init with everyone's read 
        // size. after read_verify, di.dims is populated with the full array 
        // shape. so the final set is needed to get a proper array shape for 
        // each rank
        //
        // collective scalar read is prohibited
        //
        // in non-collective vector read, all_dim0 is inited with all 0, and
        // di.dims are set to the dims of the full array. In this case, set the
        // root rank of all_dims0 to the full dim0, and leave all other ranks 
        // to 0. Later in read_impl it doesnt need to read the file again in
        // order to figure out the full array size (so when the size is 0, the
        // read should be skipped)
        //
        // in non-collective scalar read, all_dim0 already init with 1 in the
        // root rank, and 0 in others. so no need to do any more setups.
        //
        if (collective) di.dims[0] = all_dims_0[mpi_rank];
        else if(data_rank) all_dims_0[root_rank] = di.dims[0];
    }

    static void 
    read_impl( 
            Hdf5_handler const& file,
            std::string const& name,
            std::vector<hsize_t> const& all_dims_0,
            syn::data_info_t const& di,
            Commxx const& comm )
    {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        // data rank
        auto data_rank = di.dims.size();

        // offsets for each rank (offsets of dim0 in the combined array)
        std::vector<hsize_t> offsets(mpi_size, 0);
        for (int r=0; r<mpi_size-1; ++r) offsets[r+1] = offsets[r] + all_dims_0[r];

        // total dim0
        hsize_t dim0 = offsets[mpi_size-1] + all_dims_0[mpi_size-1];

#ifdef USE_PARALLEL_HDF5

        if (!file.valid())
            throw std::runtime_error("invalid file handler");

        Hdf5_handler dset   = H5Dopen(file, name.c_str(), H5P_DEFAULT);
        Hdf5_handler fspace = H5Dget_space(dset);
        Hdf5_handler mspace = H5Screate_simple(data_rank, di.dims.data(), NULL);

        if (data_rank)
        {
            // do not read if the file data is empty
            if (dim0 == 0) return;

            // select the slab
            auto offset = std::vector<hsize_t>(data_rank, 0);
            offset[0] = offsets[mpi_rank];

            auto res = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, 
                    offset.data(), NULL, di.dims.data(), NULL);
            
            if (res < 0) throw Hdf5_exception("error at select hyperslab");
        }

        // collective write
        Hdf5_handler plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // read
        auto res = H5Dread(dset, di.atomic_type, mspace, fspace, plist_id, (void*)di.ptr);
        if (res < 0) throw Hdf5_exception("error read");
        //if (res < 0) H5Eprint(H5E_DEFAULT, stdout);

#else
        // TODO: serial hdf5 impl
        // ....
#endif
    }

};

#endif
