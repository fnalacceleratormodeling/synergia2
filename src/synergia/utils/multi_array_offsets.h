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

#endif /* MULTI_ARRAY_OFFSETS_H_ */
