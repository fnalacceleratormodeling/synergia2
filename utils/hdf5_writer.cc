#include "hdf5_writer.h"

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name) :
        file(file), name(name), data_rank(1), dims(1)
    {
        atomic_type = h5_atomic_typename<double > ();
    }

template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name) :
        file(file), name(name), data_rank(2), dims(2)
    {
        atomic_type = h5_atomic_typename<double > ();
    }

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name) :
        file(file), name(name), data_rank(3), dims(3)
    {
        atomic_type = h5_atomic_typename<double > ();
    }

template<>
    void
    Hdf5_writer<MArray1d_ref >::update_dims(MArray1d_ref & data)
    {
        dims.at(0) = data.shape()[0];
    }

template<>
    void
    Hdf5_writer<MArray2d_ref >::update_dims(MArray2d_ref & data)
    {
        dims.at(0) = data.shape()[0];
        dims.at(1) = data.shape()[1];
    }
template<>
    void
    Hdf5_writer<MArray3d_ref >::update_dims(MArray3d_ref & data)
    {
        dims.at(0) = data.shape()[0];
        dims.at(1) = data.shape()[1];
        dims.at(2) = data.shape()[2];
    }

template<>
    void *
    Hdf5_writer<MArray1d_ref >::get_data_ptr(MArray1d_ref & data)
    {
        return data.origin();
    }

template<>
    void *
    Hdf5_writer<MArray2d_ref >::get_data_ptr(MArray2d_ref & data)
    {
        return data.origin();
    }
template<>
    void *
    Hdf5_writer<MArray3d_ref >::get_data_ptr(MArray3d_ref & data)
    {
        return data.origin();
    }

