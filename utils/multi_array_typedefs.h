#ifndef MULTI_ARRAY_TYPEDEFS_H_
#define MULTI_ARRAY_TYPEDEFS_H_

#include "boost/multi_array.hpp"
typedef boost::multi_array<double, 1 > MArray1d;
typedef boost::multi_array_ref<double, 1 > MArray1d_ref;
typedef boost::const_multi_array_ref<double, 1 > Const_MArray1d_ref;

typedef boost::multi_array<double, 2 > MArray2d;
typedef boost::multi_array_ref<double, 2 > MArray2d_ref;
typedef boost::const_multi_array_ref<double, 2 > Const_MArray2d_ref;

typedef boost::multi_array<double, 3 > MArray3d;
typedef boost::multi_array_ref<double, 3 > MArray3d_ref;
typedef boost::const_multi_array_ref<double, 3 > Const_MArray3d_ref;

#endif /* MULTI_ARRAY_TYPEDEFS_H_ */
