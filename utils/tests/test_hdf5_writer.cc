#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/hdf5_writer.h"
#include "utils/multi_array_typedefs.h"

// jfa: these are bad tests because they require the user
// to manually inspect the output files.

BOOST_AUTO_TEST_CASE(construct_integer)
{
    hid_t file = H5Fcreate("integer.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    Hdf5_writer<int > writer(file, "i");
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(integer)
{
    hid_t file = H5Fcreate("integer.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    Hdf5_writer<int > writer(file, "i");
    for (int i = 0; i < 5; ++i) {
        writer.append(i);
    }
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(double_float)
{
    hid_t file = H5Fcreate("double_float.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    Hdf5_writer<double > writer(file, "double_float");
    for (int i = 0; i < 5; ++i) {
        double x = 1.1 * i;
        writer.append(x);
    }
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array1d)
{
    hid_t file = H5Fcreate("array1d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim = 6;
    MArray1d a(boost::extents[dim]);
    Hdf5_writer<MArray1d_ref > writer(file, "array1d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim; ++j) {
            a[j] = 1.1 * j + 10 * i;
        }
        writer.append(a);
    }
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array2d)
{
    hid_t file = H5Fcreate("array2d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_writer<MArray2d_ref > writer(file, "array2d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                a[j][k] = 1.1 * k + 10 * j + 100 * i;
            }
        }
        writer.append(a);
    }
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array3d)
{
    hid_t file = H5Fcreate("array3d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray3d a(boost::extents[dim1][dim2][dim3]);
    Hdf5_writer<MArray3d_ref > writer(file, "array3d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                for (int l = 0; l < dim2; ++l) {
                    a[j][k][l] = 1.1 * l + 10 * k + 100 * j + 1000 * i;
                }
            }
        }
        writer.append(a);
    }
    H5Fclose(file);
}
