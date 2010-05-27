#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/hdf5_writer.h"
#include "utils/multi_array_typedefs.h"

BOOST_AUTO_TEST_CASE(integer)
{
    hid_t file = H5Fcreate("justinteger.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    Hdf5_writer<int > writer(file, "i");
    int data = 7;
    writer.write(data);
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(double_float)
{
    hid_t file = H5Fcreate("justdouble.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    Hdf5_writer<double > writer(file, "double_float");
    double data = 2.71828;
    writer.write(data);
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array1d)
{
    hid_t file = H5Fcreate("justarray1d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim = 6;
    MArray1d a(boost::extents[dim]);
    Hdf5_writer<MArray1d_ref > writer(file, "array1d");
    for (int j = 0; j < dim; ++j) {
        a[j] = 1.1 * j;
    }
    writer.write(a);
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array2d)
{
    hid_t file = H5Fcreate("justarray2d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_writer<MArray2d_ref > writer(file, "array2d");
    for (int j = 0; j < dim1; ++j) {
        for (int k = 0; k < dim2; ++k) {
            a[j][k] = 1.1 * k + 10 * j;
        }
    }
    writer.write(a);
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array3d)
{
    hid_t file = H5Fcreate("justarray3d.h5", H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray3d a(boost::extents[dim1][dim2][dim3]);
    Hdf5_writer<MArray3d_ref > writer(file, "array3d");
    for (int j = 0; j < dim1; ++j) {
        for (int k = 0; k < dim2; ++k) {
            for (int l = 0; l < dim2; ++l) {
                a[j][k][l] = 1.1 * l + 10 * k + 100 * j;
            }
        }
    }
    writer.write(a);
    H5Fclose(file);
}

