#ifndef MULTI_ARRAY_CHECK_EQUAL_H_
#define MULTI_ARRAY_CHECK_EQUAL_H_

inline void
multi_array_check_equal(Const_MArray1d_ref const& a,
        Const_MArray1d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    for (unsigned int i = 0; i < a.shape()[0]; ++i) {
        if (a[i] == 0.0) {
            BOOST_CHECK_SMALL(b[i], tolerance);
        } else {
            BOOST_CHECK_CLOSE(a[i], b[i], tolerance);
        }
    }
}

inline void
multi_array_check_equal(Const_MArray2d_ref const& a,
        Const_MArray2d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    for (unsigned int i = 0; i < a.shape()[0]; ++i) {
        for (unsigned int j = 0; j < a.shape()[1]; ++j) {
            if (a[i][j] == 0.0) {
                BOOST_CHECK_SMALL(b[i][j], tolerance);
            } else {
                BOOST_CHECK_CLOSE(a[i][j], b[i][j], tolerance);
            }
        }
    }
}

inline void
multi_array_check_equal(Const_MArray3d_ref const& a,
        Const_MArray3d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    BOOST_CHECK_EQUAL(a.shape()[2], b.shape()[2]);
    for (unsigned int i = 0; i < a.shape()[0]; ++i) {
        for (unsigned int j = 0; j < a.shape()[1]; ++j) {
            for (unsigned int k = 0; k < a.shape()[2]; ++k) {
                if (a[i][j][k] == 0.0) {
                    BOOST_CHECK_SMALL(b[i][j][k], tolerance);
                } else {
                    BOOST_CHECK_CLOSE(a[i][j][k], b[i][j][k], tolerance);
                }
            }
        }
    }
}

#endif /* MULTI_ARRAY_CHECK_EQUAL_H_ */
