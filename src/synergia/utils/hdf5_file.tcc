#ifndef HDF5_FILE_TCC_
#define HDF5_FILE_TCC_

#include "synergia/utils/hdf5_misc.h"
#include <boost/type_traits.hpp>

template<typename T>
    void
    Hdf5_file::write(T const& data, std::string const& name)
    {
        Hdf5_writer<T> (h5file.hid, name).write(data);
    }

template<typename T>
    void
    Hdf5_file::write(T const * data, size_t len, std::string const & name)
    {
        Hdf5_writer<T> (h5file.hid, name).write(data, len);
    }

template<typename T>
    T
    Hdf5_file::read(std::string const& name)
    {
        T retval;

        Hdf5_handler atomic_type = hdf5_atomic_data_type<T>();
        Hdf5_handler dataset = H5Dopen2(h5file.hid, name.c_str(), H5P_DEFAULT);

        std::vector<hsize_t > dims(1);
        dims.at(0) = 1;
        int data_rank = 0;

        Hdf5_handler dataspace = H5Screate_simple(data_rank, &dims[0], NULL);
        Hdf5_handler memspace = H5Screate_simple(data_rank, &dims[0], NULL);

        herr_t res = H5Dread(dataset.hid, atomic_type.hid, memspace.hid, 
                dataspace.hid, H5P_DEFAULT, &retval);

        if (res < 0) 
            throw Hdf5_exception("Error at reading Dataset " + name + " from HDF5 file");

        return retval;
    }

template<>
    int *
    Hdf5_file::read<int *    >(std::string const& name);
template<>
    double *
    Hdf5_file::read<double * >(std::string const& name);
template<>
    MArray1d
    Hdf5_file::read<MArray1d >(std::string const& name);
template<>
    MArray2d
    Hdf5_file::read<MArray2d >(std::string const& name);
template<>
    MArray3d
    Hdf5_file::read<MArray3d >(std::string const& name);
template<>
    MArray1i
    Hdf5_file::read<MArray1i >(std::string const& name);

#endif /* HDF5_FILE_TCC_ */
