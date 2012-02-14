#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_check_equal.h"

BOOST_AUTO_TEST_CASE(construct)
{
    Hdf5_file file("hdf5_file_empty.h5", Hdf5_file::truncate);
}

BOOST_AUTO_TEST_CASE(write_data)
{
    Hdf5_file file("hdf5_file_write.h5", Hdf5_file::truncate);
    int int_data = 7;
    file.write(int_data, "int_data");
    double double_data = 2.71828;
    file.write(double_data, "double_data");
    const int dim = 6;
    MArray1d a1d(boost::extents[dim]);
    for (int j = 0; j < dim; ++j) {
        a1d[j] = 1.1 * j;
    }
    file.write(a1d, "a1d");
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray3d a3d(boost::extents[dim1][dim2][dim3]);
    for (int j = 0; j < dim1; ++j) {
        for (int k = 0; k < dim2; ++k) {
            for (int l = 0; l < dim2; ++l) {
                a3d[j][k][l] = 1.1 * l + 10 * k + 100 * j;
            }
        }
    }
    file.write(a3d, "a3d");
}

const double tolerance = 1.0e-13;
BOOST_AUTO_TEST_CASE(read_write_data)
{
    const char * filename = "hdf5_file_read_write.h5";
    const char * int_label = "int_data";
    const char * double_label = "double_data";
    const char * array1d_label = "array1d_data";
    const char * array2d_label = "array2d_data";
    const char * array3d_label = "array3d_data";
    int int_data = 7;
    double double_data = 2.71828;
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray1d a1d(boost::extents[dim1]);
    MArray2d a2d(boost::extents[dim1][dim2]);
    MArray3d a3d(boost::extents[dim1][dim2][dim3]);
    for (int j = 0; j < dim1; ++j) {
        a1d[j] = 100 * j;
        for (int k = 0; k < dim2; ++k) {
            a2d[j][k] = 10 * k + 100 * j;
            for (int l = 0; l < dim2; ++l) {
                a3d[j][k][l] = 1.1 * l + 10 * k + 100 * j;
            }
        }
    }

    {
        Hdf5_file write_file(filename, Hdf5_file::truncate);
        write_file.write(int_data, int_label);
        write_file.write(double_data, double_label);
        write_file.write(a1d, array1d_label);
        write_file.write(a2d, array2d_label);
        write_file.write(a3d, array3d_label);
    }

    {
        Hdf5_file read_file(filename, Hdf5_file::read_only);
        int int_read = read_file.read<int > (int_label);
        BOOST_CHECK_EQUAL(int_read, int_data);
        double double_read = read_file.read<double > (double_label);
        BOOST_CHECK_CLOSE(double_read, double_data, tolerance);
        MArray1d a1d_read(read_file.read<MArray1d > (array1d_label));
        multi_array_check_equal(a1d_read, a1d, tolerance);
        MArray2d a2d_read(read_file.read<MArray2d > (array2d_label));
        multi_array_check_equal(a2d_read, a2d, tolerance);
        MArray3d a3d_read(read_file.read<MArray3d > (array3d_label));
        multi_array_check_equal(a3d_read, a3d, tolerance);
    }
}

