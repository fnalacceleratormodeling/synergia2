#ifndef HDF5_MISC_H_
#define HDF5_MISC_H_

#include <stdexcept>
#include "H5Cpp.h"

// The generic (T) version of h5_atomic_data_type is undefined.
// Only versions with specializations will compile.
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

#endif /* HDF5_MISC_H_ */
