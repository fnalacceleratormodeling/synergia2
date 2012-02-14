#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/hdf5_serial_writer.h"
#include "synergia/utils/multi_array_typedefs.h"

// jfa: these are bad tests because they require the user
// to manually inspect the output files.
BOOST_AUTO_TEST_CASE(construct_integer)
{
    Hdf5_file_sptr file_sptr(new Hdf5_file("integer.h5", Hdf5_file::truncate));
    Hdf5_serial_writer<int > writer(file_sptr, "i");
}

BOOST_AUTO_TEST_CASE(integer)
{
    Hdf5_file_sptr file_sptr(new Hdf5_file("integer.h5", Hdf5_file::truncate));
    Hdf5_serial_writer<int > writer(file_sptr, "i");
    for (int i = 0; i < 5; ++i) {
        writer.append(i);
    }
}

BOOST_AUTO_TEST_CASE(double_float)
{
    Hdf5_file_sptr file_sptr(
            new Hdf5_file("double_float.h5", Hdf5_file::truncate));
    Hdf5_serial_writer<double > writer(file_sptr, "double_float");
    for (int i = 0; i < 5; ++i) {
        double x = 1.1 * i;
        writer.append(x);
    }
}

BOOST_AUTO_TEST_CASE(array1d)
{
    Hdf5_file_sptr file_sptr(new Hdf5_file("array1d.h5", Hdf5_file::truncate));
    const int dim = 6;
    MArray1d a(boost::extents[dim]);
    Hdf5_serial_writer<MArray1d_ref > writer(file_sptr, "array1d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim; ++j) {
            a[j] = 1.1 * j + 10 * i;
        }
        writer.append(a);
    }
}

BOOST_AUTO_TEST_CASE(array2d)
{
    Hdf5_file_sptr file_sptr(new Hdf5_file("array2d.h5", Hdf5_file::truncate));
    const int dim1 = 2;
    const int dim2 = 3;
    MArray2d a(boost::extents[dim1][dim2]);
    Hdf5_serial_writer<MArray2d_ref > writer(file_sptr, "array2d");
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < dim1; ++j) {
            for (int k = 0; k < dim2; ++k) {
                a[j][k] = 1.1 * k + 10 * j + 100 * i;
            }
        }
        writer.append(a);
    }
}

BOOST_AUTO_TEST_CASE(array3d)
{
    Hdf5_file_sptr file_sptr(new Hdf5_file("array3d.h5", Hdf5_file::truncate));
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray3d a(boost::extents[dim1][dim2][dim3]);
    Hdf5_serial_writer<MArray3d_ref > writer(file_sptr, "array3d");
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
}
