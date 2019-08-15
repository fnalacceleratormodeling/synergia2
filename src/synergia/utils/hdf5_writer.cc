#include "hdf5_writer.h"

namespace 
{
    template <typename T>
    void storage_order_writer_helper(hid_t file_ptr, std::string const& name, T const& data)
    {
        int storage_order;

        if (data.storage_order() == boost::c_storage_order()) 
        {
            storage_order = Hdf5_writer<T>::c_storage_order;
        } 
        else if (data.storage_order() == boost::fortran_storage_order()) 
        {
            storage_order = Hdf5_writer<T>::fortran_storage_order;
        } 
        else 
        {
            throw(std::runtime_error("Hdf5_writer: unsupported storage order"));
        }

        Hdf5_writer<int> storage_order_writer(file_ptr, name + "_storage_order");
        storage_order_writer.write(storage_order);
    }
}

template<>
    void
    Hdf5_writer<MArray2d >::write_storage_order(MArray2d const& data)
    {
        storage_order_writer_helper(file_ptr, name, data);
    }

template<>
    void
    Hdf5_writer<MArray2d_ref >::write_storage_order(MArray2d_ref const& data)
    {
        storage_order_writer_helper(file_ptr, name, data);
    }

    template<>
    void
    Hdf5_writer<MArray3d >::write_storage_order(MArray3d const& data)
    {
        storage_order_writer_helper(file_ptr, name, data);
    }

    template<>
    void
    Hdf5_writer<MArray3d_ref >::write_storage_order(MArray3d_ref const& data)
    {
        storage_order_writer_helper(file_ptr, name, data);
    }

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(hid_t file_ptr,
            std::string const& name) :
        data_rank(1), dims(1), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<double>())
    {
    }

template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(hid_t file_ptr,
            std::string const& name) :
        data_rank(2), dims(2), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<double>())
    {
    }

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(hid_t file_ptr,
            std::string const& name) :
        data_rank(3), dims(3), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<double>())
    {
    }

template<>
    void
    Hdf5_writer<MArray1d_ref >::update_dims(MArray1d_ref const& data)
    {
        dims.at(0) = data.shape()[0];
    }

template<>
    void
    Hdf5_writer<MArray2d_ref >::update_dims(MArray2d_ref const& data)
    {
        dims.at(0) = data.shape()[0];
        dims.at(1) = data.shape()[1];
    }
template<>
    void
    Hdf5_writer<MArray3d_ref >::update_dims(MArray3d_ref const& data)
    {
        dims.at(0) = data.shape()[0];
        dims.at(1) = data.shape()[1];
        dims.at(2) = data.shape()[2];
    }

template<>
    const void *
    Hdf5_writer<MArray1d_ref >::get_data_ptr(MArray1d_ref const& data)
    {
        return data.origin();
    }

template<>
    const void *
    Hdf5_writer<MArray2d_ref >::get_data_ptr(MArray2d_ref const& data)
    {
        return data.origin();
    }
template<>
    const void *
    Hdf5_writer<MArray3d_ref >::get_data_ptr(MArray3d_ref const& data)
    {
        return data.origin();
    }

template<>
    Hdf5_writer<MArray1d >::Hdf5_writer(hid_t file_ptr,
            std::string const& name) :
        data_rank(1), dims(1), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<double>())
    {
    }

template<>
    Hdf5_writer<MArray2d >::Hdf5_writer(hid_t file_ptr,
            std::string const& name) :
        data_rank(2), dims(2), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<double>())
    {
    }

template<>
    Hdf5_writer<MArray3d >::Hdf5_writer(hid_t file_ptr,
            std::string const& name) :
        data_rank(3), dims(3), name(name), file_ptr(file_ptr),
        atomic_type(hdf5_atomic_data_type<double>())
    {
    }

template<>
    void
    Hdf5_writer<MArray1d >::update_dims(MArray1d const& data)
    {
        dims.at(0) = data.shape()[0];
    }

template<>
    void
    Hdf5_writer<MArray2d >::update_dims(MArray2d const& data)
    {
        dims.at(0) = data.shape()[0];
        dims.at(1) = data.shape()[1];
    }
template<>
    void
    Hdf5_writer<MArray3d >::update_dims(MArray3d const& data)
    {
        dims.at(0) = data.shape()[0];
        dims.at(1) = data.shape()[1];
        dims.at(2) = data.shape()[2];
    }

template<>
    const void *
    Hdf5_writer<MArray1d >::get_data_ptr(MArray1d const& data)
    {
        return data.origin();
    }

template<>
    const void *
    Hdf5_writer<MArray2d >::get_data_ptr(MArray2d const& data)
    {
        return data.origin();
    }
template<>
    const void *
    Hdf5_writer<MArray3d >::get_data_ptr(MArray3d const& data)
    {
        return data.origin();
    }

