#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/hdf5_serial_writer.h"
#include "synergia/utils/multi_array_typedefs.h"

// jfa: these are bad tests because they require the user
// to manually inspect the output files.

BOOST_AUTO_TEST_CASE(construct_integer)
{
    H5::H5File file("integer.h5", H5F_ACC_TRUNC);
    Hdf5_serial_writer<int > writer(file, "i");
    file.close();
}

BOOST_AUTO_TEST_CASE(integer)
{
    H5::H5File file("integer.h5", H5F_ACC_TRUNC);
    Hdf5_serial_writer<int > writer(file, "i");
    for (int i = 0; i < 5; ++i) {
        writer.append(i);
    }
    file.close();
}

BOOST_AUTO_TEST_CASE(double_float)
{
    H5::H5File file("double_float.h5", H5F_ACC_TRUNC);
    Hdf5_serial_writer<double > writer(file, "double_float");
    for (int i = 0; i < 5; ++i) {
        double x = 1.1 * i;
        writer.append(x);
    }
    file.close();
}

BOOST_AUTO_TEST_CASE(array1d)
{
    H5::H5File file("array1d.h5", H5F_ACC_TRUNC);
    const int dim = 6;
    MArray1d a(boost::extents[dim]);
    Hdf5_serial_writer<MArray1d_ref > writer(file, "array1d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim; ++j) {
            a[j] = 1.1 * j + 10 * i;
        }
        writer.append(a);
    }
    file.close();
}

BOOST_AUTO_TEST_CASE(array2d)
{
    H5::H5File file("array2d.h5", H5F_ACC_TRUNC);
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_serial_writer<MArray2d_ref > writer(file, "array2d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                a[j][k] = 1.1 * k + 10 * j + 100 * i;
            }
        }
        writer.append(a);
    }
    file.close();
}

BOOST_AUTO_TEST_CASE(array3d)
{
    H5::H5File file("array3d.h5", H5F_ACC_TRUNC);
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray3d a(boost::extents[dim1][dim2][dim3]);
    Hdf5_serial_writer<MArray3d_ref > writer(file, "array3d");
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
    file.close();
}
