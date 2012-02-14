#ifndef HDF5_FILE_TCC_
#define HDF5_FILE_TCC_

#include "synergia/utils/hdf5_misc.h"

template<typename T>
    void
    Hdf5_file::write(T const& data, std::string const& name)
    {
        Hdf5_writer<T > (h5file_ptr, name).write(data);
    }

template<typename T>
    T
    Hdf5_file::read(std::string const& name)
    {
        T retval;
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<T > ();
        std::vector<hsize_t > dims(1);
        dims.at(0) = 1;
        int data_rank = 0;
        DataSpace dataspace(data_rank, &dims[0]);
        DataSpace memspace(data_rank, &dims[0]);
        T * data_out = &retval;
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
    }

template<>
    MArray1d
    Hdf5_file::read<MArray1d >(std::string const& name);
template<>
    MArray2d
    Hdf5_file::read<MArray2d >(std::string const& name);
template<>
    MArray3d
    Hdf5_file::read<MArray3d >(std::string const& name);

#endif /* HDF5_FILE_TCC_ */
