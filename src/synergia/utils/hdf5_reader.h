#ifndef HDF5_READER_H
#define HDF5_READER_H

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/multi_array_typedefs.h"

class Hdf5_reader
{
private:

    Hdf5_reader() { };


public:

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
        return T();
    }

    template<class T>
    static T
    read( Hdf5_handler const& file,
          std::string const& name,
          size_t len,
          Commxx const& comm,
          int root_rank )
    {
        return T();
    }

    // read 1d array
    template<class T>
    static void
    read( Hdf5_handler const& file,
          std::string const& name,
          T const* data, 
          Commxx const& comm,
          int root_rank )
    {
    }
 
    // collective read a 1d array
    template<class T>
    static void
    read( Hdf5_handler const& file,
          std::string const& name,
          T const* data, 
          size_t len,
          Commxx const& comm,
          int root_rank )
    {
#if 0
        auto data_ptr = data;
        auto data_dims = std::vector<hsize_t>({len});
        auto atomic = hdf5_atomic_data_type<T>();

        auto all_dim0 = collective_dims();

        read_impl(file, name, data_ptr, data_dims, all_dim0, atomic, comm);
#endif
    }

    void read_impl()
    {
#if 0
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
#else
#endif

#endif
    }

};

#endif
