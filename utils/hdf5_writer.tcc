#include "utils/multi_array_typedefs.h"
#include "utils/hdf5_utils.h"
#include <stdexcept>
#include <iostream>

template<typename T>
    Hdf5_writer<T >::Hdf5_writer(hid_t & file, std::string const& name) :
        file(file), name(name), data_rank(0), dims(1)
    {
        atomic_type = h5_atomic_typename<T > ();
    }

template<typename T>
    void
    Hdf5_writer<T >::update_dims(T & data)
    {
        dims.at(0) = 1;
    }

template<typename T>
    void *
    Hdf5_writer<T >::get_data_ptr(T & data)
    {
        return &data;
    }

template<typename T>
    void
    Hdf5_writer<T >::write(T & data)
    {
        update_dims(data);
        hid_t dataspace_id = H5Screate_simple(data_rank, &dims[0], NULL);
        hid_t dataset_id = H5Dcreate(file, name.c_str(), atomic_type,
                dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t status = H5Dwrite(dataset_id, atomic_type, H5S_ALL, H5S_ALL,
                H5P_DEFAULT, get_data_ptr(data));
        h5_error_check(status);
        status = H5Dclose(dataset_id);
        h5_error_check(status);
        status = H5Sclose(dataspace_id);
        h5_error_check(status);
    }

template<typename T>
    Hdf5_writer<T >::~Hdf5_writer()
    {
    }

#include "hdf5_writer.h"

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name);

template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name);

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name);

template<>
    void
    Hdf5_writer<MArray1d_ref >::update_dims(MArray1d_ref & data);

template<>
    void
    Hdf5_writer<MArray2d_ref >::update_dims(MArray2d_ref & data);

template<>
    void
    Hdf5_writer<MArray3d_ref >::update_dims(MArray3d_ref & data);

template<>
    void *
    Hdf5_writer<MArray1d_ref >::get_data_ptr(MArray1d_ref & data);

template<>
    void *
    Hdf5_writer<MArray2d_ref >::get_data_ptr(MArray2d_ref & data);

template<>
    void *
    Hdf5_writer<MArray3d_ref >::get_data_ptr(MArray3d_ref & data);

