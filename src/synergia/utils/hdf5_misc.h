#ifndef HDF5_MISC_H_
#define HDF5_MISC_H_

#include <stdexcept>
#include "hdf5.h"

inline void
hdf5_error_check(hid_t status, const char * message = "hdf5 error")
{
    if (status != 0) {
        throw(std::runtime_error(message));
    }
}

// The generic (T) version of h5_atomic_typenameis undefined.
// Only versions with specializations will compile.
template<typename T>
    inline hid_t
    hdf5_atomic_typename();

template<>
    inline hid_t
    hdf5_atomic_typename<int > ()
    {
        return H5T_NATIVE_INT;
    }

template<>
    inline hid_t
    hdf5_atomic_typename<double > ()
    {
        return H5T_NATIVE_DOUBLE;
    }

#endif /* HDF5_MISC_H_ */
