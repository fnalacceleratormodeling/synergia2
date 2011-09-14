#ifndef MULTI_ARRAY_OFFSETS_H_
#define MULTI_ARRAY_OFFSETS_H_
#include "synergia/utils/multi_array_typedefs.h"

inline
double *
multi_array_offset(MArray3d_ref & array, int i, int j, int k)
{
    return array.origin() + i * array.strides()[0] + j * array.strides()[1] + k
            * array.strides()[2];
}

inline
const double *
multi_array_offset(MArray3d_ref const& array, int i, int j, int k)
{
    return array.origin() + i * array.strides()[0] + j * array.strides()[1] + k
            * array.strides()[2];
}

inline
std::complex<double > *
multi_array_offset(MArray3dc_ref & array, int i, int j, int k)
{
    return array.origin() + i * array.strides()[0] + j * array.strides()[1] + k
            * array.strides()[2];
}

inline
double *
multi_array_offset(MArray2d_ref & array, int i, int j)
{
    return array.origin() + i * array.strides()[0] + j * array.strides()[1];
}

inline
const double *
multi_array_offset(MArray2d_ref const& array, int i, int j)
{
    return array.origin() + i * array.strides()[0] + j * array.strides()[1];
}

inline
std::complex<double > *
multi_array_offset(MArray2dc_ref & array, int i, int j)
{
    return array.origin() + i * array.strides()[0] + j * array.strides()[1];
}

inline
double *
multi_array_offset(MArray1d_ref & array, int i)
{
    return array.origin() + i * array.strides()[0];
}

inline
const double *
multi_array_offset(MArray1d_ref const& array, int i)
{
    return array.origin() + i * array.strides()[0];
}

inline
std::complex<double > *
multi_array_offset(MArray1dc_ref & array, int i)
{
    return array.origin() + i * array.strides()[0];
}

#endif /* MULTI_ARRAY_OFFSETS_H_ */
