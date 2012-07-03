#ifndef MULTI_ARRAY_CHECK_EQUAL_H_
#define MULTI_ARRAY_CHECK_EQUAL_H_
#include "synergia/utils/floating_point.h"

inline void
multi_array_check_equal(Const_MArray1d_ref const& a,
        Const_MArray1d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    for (unsigned int i = a.index_bases()[0]; i < a.index_bases()[0] + a.shape()[0]; ++i) {
        BOOST_CHECK_MESSAGE(floating_point_equal(a[i], b[i], tolerance), "a["
                << i << "] = " << a[i] << ", b[" << i << "] = " << b[i]
                << ", a-b = " << a[i] - b[i] << ", tolerance = " << tolerance);
    }
}

inline void
multi_array_check_equal(Const_MArray2d_ref const& a,
        Const_MArray2d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.index_bases()[1], b.index_bases()[1]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    for (unsigned int i = a.index_bases()[0]; i < a.index_bases()[0] + a.shape()[0]; ++i) {
        for (unsigned int j = a.index_bases()[1]; j < a.index_bases()[1] + a.shape()[1]; ++j) {
            BOOST_CHECK_MESSAGE(floating_point_equal(a[i][j], b[i][j],
                    tolerance), "a[" << i << "][" << j << "] = " << a[i][j]
                    << ", b[" << i << "][" << j << "] = " << b[i][j]
                    << ", a-b = " << a[i][j] - b[i][j] << ", tolerance = "
                    << tolerance);
        }
    }
}

inline void
multi_complex_array_check_equal(Const_MArray2dc_ref const& a,
        Const_MArray2dc_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.index_bases()[1], b.index_bases()[1]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    for (unsigned int i = a.index_bases()[0];
            i < a.index_bases()[0] + a.shape()[0]; ++i) {
        for (unsigned int j = a.index_bases()[1];
                j < a.index_bases()[1] + a.shape()[1]; ++j) {
            //if (abs(a[i][j]) < tolerance) {
            BOOST_CHECK_MESSAGE(
                    complex_floating_point_equal(a[i][j], b[i][j],
                    tolerance), "a[" << i << "][" << j << "] = " << a[i][j]
                    << ", b[" << i << "][" << j << "] = " << b[i][j]
                    << ", |a-b| = " << abs(a[i][j] - b[i][j])
                    << ", tolerance = " << tolerance);
            //} else {
            //BOOST_CHECK_MESSAGE(complex_floating_point_equal(a[i][j], b[i][j],
            //        tolerance), "a[" << i << "][" << j << "] = " << a[i][j]
            //        << ", b[" << i << "][" << j << "] = " << b[i][j]
            //        << ", |(a-b)/a| = " << abs((a[i][j] - b[i][j]) / a[i][j]) 
            //        << ", tolerance = " << tolerance);
            //}
        }
    }
}

inline void
multi_array_check_equal(Const_MArray3d_ref const& a,
        Const_MArray3d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.index_bases()[1], b.index_bases()[1]);
    BOOST_CHECK_EQUAL(a.index_bases()[2], b.index_bases()[2]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    BOOST_CHECK_EQUAL(a.shape()[2], b.shape()[2]);
    for (unsigned int i = a.index_bases()[0]; i < a.index_bases()[0] + a.shape()[0]; ++i) {
        for (unsigned int j = a.index_bases()[1]; j < a.index_bases()[1] + a.shape()[1]; ++j) {
            for (unsigned int k = a.index_bases()[2]; k < a.index_bases()[2]
                    + a.shape()[2]; ++k) {
                BOOST_CHECK_MESSAGE(floating_point_equal(a[i][j][k],
                        b[i][j][k], tolerance), "a[" << i << "][" << j << "]["
                        << k << "] = " << a[i][j][k] << ", b[" << i << "]["
                        << j << "][" << k << "] = " << b[i][j][k] << ", a-b = "
                        << a[i][j][k] - b[i][j][k] << ", tolerance = "
                        << tolerance);
            }
        }
    }
}

#endif /* MULTI_ARRAY_CHECK_EQUAL_H_ */
