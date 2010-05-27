#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/hdf5_chunked_array2d_writer.h"

// jfa: these are bad tests because they require the user
// to manually inspect the output files.

BOOST_AUTO_TEST_CASE(construct)
{
    hid_t file = H5Fcreate("chunkedarray2d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(file, "array2d", a);
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(manual_close)
{
    hid_t file = H5Fcreate("chunkedarray2d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(file, "array2d", a);
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(typical)
{
    hid_t file = H5Fcreate("chunkedarray2d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(file, "array2d", a);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                a[j][k] = 1.1 * k + 10 * j + 100 * (i + 1);
            }
        }
        writer.write_chunk(a);
    }
    writer.close();
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(write_to_closed)
{
    hid_t file = H5Fcreate("chunkedarray2d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_chunked_array2d_writer writer(file, "array2d", a);
    writer.close();
    bool caught_error = false;
    try {
        writer.write_chunk(a);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }
    H5Fclose(file);
    BOOST_CHECK(caught_error == true);
}

