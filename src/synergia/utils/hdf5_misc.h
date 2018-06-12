#ifndef HDF5_MISC_H_
#define HDF5_MISC_H_

#include <stdexcept>
//#include "H5Cpp.h"
#include "hdf5.h"

// The generic (T) version of h5_atomic_data_type is undefined.
// Only versions with specializations will compile.
#if 0
template<typename T>
    inline H5::DataType
    hdf5_atomic_data_type();

template<>
    inline H5::DataType
    hdf5_atomic_data_type<int > ()
    {
        return H5::PredType::NATIVE_INT;
    }

template<>
    inline H5::DataType
    hdf5_atomic_data_type<double > ()
    {
        return H5::PredType::NATIVE_DOUBLE;
    }
#endif


template<typename T>
    inline hid_t
    hdf5_atomic_data_type();

template<>
    inline hid_t
    hdf5_atomic_data_type<int > ()
    {
        return H5Tcopy(H5T_NATIVE_INT);
    }

template<>
    inline hid_t
    hdf5_atomic_data_type<double > ()
    {
        return H5Tcopy(H5T_NATIVE_DOUBLE);
    }


#endif /* HDF5_MISC_H_ */
