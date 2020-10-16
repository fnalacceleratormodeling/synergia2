#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"
#include <stdexcept>
#include <iostream>

template<typename T>
    Hdf5_writer<T >::Hdf5_writer(hid_t file_ptr, std::string const& name) :
        data_rank(0), dims(1), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<T>())
    {
    }

template<typename T>
    void
    Hdf5_writer<T >::write_storage_order(T const& data)
    {
    }

template<>
    void
    Hdf5_writer<MArray2d >::write_storage_order(MArray2d const& data);

template<>
    void
    Hdf5_writer<MArray2d_ref >::write_storage_order(MArray2d_ref const& data);

template<>
    void
    Hdf5_writer<MArray3d >::write_storage_order(MArray3d const& data);

template<>
    void
    Hdf5_writer<MArray3d_ref >::write_storage_order(MArray3d_ref const& data);

template<typename T>
    void
    Hdf5_writer<T >::update_dims(T const& data)
    {
        dims.at(0) = 1;
    }

template<typename T>
    const void *
    Hdf5_writer<T >::get_data_ptr(T const& data)
    {
        return &data;
    }

template<typename T>
    void
    Hdf5_writer<T >::write(T const& data)
    {
        write_storage_order(data);
        update_dims(data);

        Hdf5_handler dataspace = H5Screate_simple(data_rank, &dims[0], NULL);
        Hdf5_handler dataset   = H5Dcreate(file_ptr, name.c_str(), atomic_type,
                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        herr_t res = H5Dwrite(dataset, atomic_type, dataspace, H5S_ALL, H5P_DEFAULT, 
                get_data_ptr(data));

        if (res < 0) throw Hdf5_exception();
    }

template<typename T>
    void
    Hdf5_writer<T >::write(T const * data, size_t len)
    {
        dims.at(0) = len;
        data_rank = 1;

        Hdf5_handler dataspace = H5Screate_simple(data_rank, &dims[0], NULL);
        Hdf5_handler dataset   = H5Dcreate(file_ptr, name.c_str(), atomic_type,
                dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        herr_t res = H5Dwrite(dataset, atomic_type, dataspace, H5S_ALL, H5P_DEFAULT, 
                (const void *)data);

        if (res < 0) throw Hdf5_exception();
    }

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(hid_t file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(hid_t file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(hid_t file_ptr,
            std::string const& name);

template<>
    void
    Hdf5_writer<MArray1d_ref >::update_dims(MArray1d_ref const& data);

template<>
    void
    Hdf5_writer<MArray2d_ref >::update_dims(MArray2d_ref const& data);

template<>
    void
    Hdf5_writer<MArray3d_ref >::update_dims(MArray3d_ref const& data);

template<>
    const void *
    Hdf5_writer<MArray1d_ref >::get_data_ptr(MArray1d_ref const& data);

template<>
    const void *
    Hdf5_writer<MArray2d_ref >::get_data_ptr(MArray2d_ref const& data);

template<>
    const void *
    Hdf5_writer<MArray3d_ref >::get_data_ptr(MArray3d_ref const& data);

template<>
    Hdf5_writer<MArray1d >::Hdf5_writer(hid_t file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray2d >::Hdf5_writer(hid_t file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray3d >::Hdf5_writer(hid_t file_ptr,
            std::string const& name);

template<>
    void
    Hdf5_writer<MArray1d >::update_dims(MArray1d const& data);

template<>
    void
    Hdf5_writer<MArray2d >::update_dims(MArray2d const& data);

template<>
    void
    Hdf5_writer<MArray3d >::update_dims(MArray3d const& data);

template<>
    const void *
    Hdf5_writer<MArray1d >::get_data_ptr(MArray1d const& data);

template<>
    const void *
    Hdf5_writer<MArray2d >::get_data_ptr(MArray2d const& data);

template<>
    const void *
    Hdf5_writer<MArray3d >::get_data_ptr(MArray3d const& data);
