#ifndef MULTI_ARRAY_CHECK_EQUAL_H_
#define MULTI_ARRAY_CHECK_EQUAL_H_

inline void
multi_array_check_equal(Const_MArray1d_ref const& a,
        Const_MArray1d_ref const& b, double tolerance)
{
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    for (int i = a.index_bases()[0]; i < a.index_bases()[0] + a.shape()[0]; ++i) {
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
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.index_bases()[1], b.index_bases()[1]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    for (int i = a.index_bases()[0]; i < a.index_bases()[0] + a.shape()[0]; ++i) {
        for (int j = a.index_bases()[1]; j < a.index_bases()[1] + a.shape()[1]; ++j) {
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
    BOOST_CHECK_EQUAL(a.index_bases()[0], b.index_bases()[0]);
    BOOST_CHECK_EQUAL(a.index_bases()[1], b.index_bases()[1]);
    BOOST_CHECK_EQUAL(a.index_bases()[2], b.index_bases()[2]);
    BOOST_CHECK_EQUAL(a.shape()[0], b.shape()[0]);
    BOOST_CHECK_EQUAL(a.shape()[1], b.shape()[1]);
    BOOST_CHECK_EQUAL(a.shape()[2], b.shape()[2]);
    for (int i = a.index_bases()[0]; i < a.index_bases()[0] + a.shape()[0]; ++i) {
        for (int j = a.index_bases()[1]; j < a.index_bases()[1] + a.shape()[1]; ++j) {
            for (int k = a.index_bases()[2]; k < a.index_bases()[2]
                    + a.shape()[2]; ++k) {
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
