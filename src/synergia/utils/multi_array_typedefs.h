#ifndef MULTI_ARRAY_TYPEDEFS_H_
#define MULTI_ARRAY_TYPEDEFS_H_

#include <complex>
#include "boost/multi_array.hpp"

typedef boost::multi_array<double, 1 > MArray1d;
typedef boost::multi_array_ref<double, 1 > MArray1d_ref;
typedef boost::const_multi_array_ref<double, 1 > Const_MArray1d_ref;
typedef boost::detail::multi_array::multi_array_view<double, 1 > MArray1d_view;

typedef boost::multi_array<double, 2 > MArray2d;
typedef boost::multi_array_ref<double, 2 > MArray2d_ref;
typedef boost::const_multi_array_ref<double, 2 > Const_MArray2d_ref;
typedef boost::detail::multi_array::multi_array_view<double, 2 > MArray2d_view;
typedef boost::detail::multi_array::multi_array_view<double, 2 > MArray2d_view;
typedef boost::detail::multi_array::const_multi_array_view<double, 2 >
        Const_MArray2d_view;

typedef boost::multi_array<double, 3 > MArray3d;
typedef boost::multi_array_ref<double, 3 > MArray3d_ref;
typedef boost::const_multi_array_ref<double, 3 > Const_MArray3d_ref;
typedef boost::detail::multi_array::multi_array_view<double, 3 > MArray3d_view;

typedef boost::multi_array<std::complex<double >, 1 > MArray1dc;
typedef boost::multi_array_ref<std::complex<double >, 1 > MArray1dc_ref;
typedef boost::const_multi_array_ref<std::complex<double >, 1 >
        Const_MArray1dc_ref;
typedef boost::detail::multi_array::multi_array_view<std::complex<double >, 1 >
        MArray1dc_view;

typedef boost::multi_array<std::complex<double >, 2 > MArray2dc;
typedef boost::multi_array_ref<std::complex<double >, 2 > MArray2dc_ref;
typedef boost::const_multi_array_ref<std::complex<double >, 2 >
        Const_MArray2dc_ref;
typedef boost::detail::multi_array::multi_array_view<std::complex<double >, 2 >
        MArray2dc_view;

typedef boost::multi_array<std::complex<double >, 3 > MArray3dc;
typedef boost::multi_array_ref<std::complex<double >, 3 > MArray3dc_ref;
typedef boost::const_multi_array_ref<std::complex<double >, 3 >
        Const_MArray3dc_ref;
typedef boost::detail::multi_array::multi_array_view<std::complex<double >, 3 >
        MArray3dc_view;

typedef boost::multi_array_types::index_range range;
typedef boost::multi_array_types::extent_range extent_range;

#endif /* MULTI_ARRAY_TYPEDEFS_H_ */
