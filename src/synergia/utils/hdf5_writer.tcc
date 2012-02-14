#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"
#include <stdexcept>
#include <iostream>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

template<typename T>
    Hdf5_writer<T >::Hdf5_writer(H5File * file_ptr, std::string const& name) :
        data_rank(0), dims(1), name(name), file_ptr(file_ptr)
    {
        atomic_type = hdf5_atomic_data_type<T > ();
    }

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
        update_dims(data);
        DataSpace dataspace(data_rank, &dims[0]);
        DataSet dataset = file_ptr->createDataSet(name.c_str(), atomic_type,
                dataspace);
        dataset.write(get_data_ptr(data), atomic_type, dataspace);

    }

template<typename T>
    Hdf5_writer<T >::~Hdf5_writer()
    {
    }

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(H5File * file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(H5File * file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(H5File * file_ptr,
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
    Hdf5_writer<MArray1d >::Hdf5_writer(H5File * file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray2d >::Hdf5_writer(H5File * file_ptr,
            std::string const& name);

template<>
    Hdf5_writer<MArray3d >::Hdf5_writer(H5File * file_ptr,
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

