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
#if 0
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
#endif
        T retval;

        hid_t dataset = H5Dopen2(h5file_ptr, name.c_str(), H5P_DEFAULT);
        hid_t atomic_type = hdf5_atomic_data_type<T>();

        std::vector<hsize_t > dims(1);
        dims.at(0) = 1;
        int data_rank = 0;

        hid_t dataspace = H5Screate_simple(data_rank, &dims[0], NULL);
        hid_t  memspace = H5Screate_simple(data_rank, &dims[0], NULL);

        H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, &retval);
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
template<>
    MArray1i
    Hdf5_file::read<MArray1i >(std::string const& name);

#endif /* HDF5_FILE_TCC_ */
