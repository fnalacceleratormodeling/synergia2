#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/hdf5_chunked_array2d_writer.h"

// jfa: these are bad tests because they require the user
// to manually inspect the output files.

BOOST_AUTO_TEST_CASE(construct)
{
    H5::H5File file("chunkedarray2d.h5", H5F_ACC_TRUNC);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(&file, "array2d", a);
    file.close();
}

BOOST_AUTO_TEST_CASE(typical)
{
    H5::H5File file("chunkedarray2d.h5", H5F_ACC_TRUNC);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(&file, "array2d", a);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                a[j][k] = 1.1 * k + 10 * j + 100 * (i + 1);
            }
        }
        writer.write_chunk(a);
    }
    file.close();
}

BOOST_AUTO_TEST_CASE(different_sizes)
{
    H5::H5File file("chunkedarray2d.h5", H5F_ACC_TRUNC);
    const int dim1 = 33;
    const int dim2 = 7;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(&file, "array2d", a);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                a[j][k] = 1.1 * k + 10 * j + 100 * (i + 1);
            }
        }
        writer.write_chunk(a);
    }
    const int dim1b = 30;
    MArray2d b(boost::extents[dim1b][dim2]);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1b; ++j) {
            for (int k = 0; k < dim2; ++k) {
                b[j][k] = 1.1 * k + 10 * j + 100 * (i + 1);
            }
        }
        writer.write_chunk(b);
    }
    file.close();
}

// ugh. no worky.
//BOOST_AUTO_TEST_CASE(typical_view)
//{
//    H5::H5File file("chunkedarray2dview.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
//            H5P_DEFAULT);
//    const int dim1 = 2;
//    const int dim2 = 3;
//    MArray2d a(boost::extents[dim1][dim2]);
//    Hdf5_chunked_array2d_writer writer(&file, "array2d", a);
//    MArray2d_view a_view = a[boost::indices[range()][range(0,2)]];
//    for (int i = 0; i < 5; ++i) {
//        for (int j = 0; j < dim1; ++j) {
//            for (int k = 0; k < dim2; ++k) {
//                a[j][k] = 1.1 * k + 10 * j + 100 * (i + 1);
//            }
//        }
//        writer.write_chunk(a_view);
//    }
//    writer.close();
//    file.close();
//}
