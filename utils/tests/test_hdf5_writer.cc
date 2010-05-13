#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/hdf5_writer.h"
#include "utils/multi_array_typedefs.h"

BOOST_AUTO_TEST_CASE(integer)
{
    hid_t file = H5Fcreate("integer.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    Hdf5_writer<int > writer(file, "i");
    for(int i=0; i<5; ++i) {
        writer.append(i);
    }
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(double_float)
{
    hid_t file = H5Fcreate("double_float.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    Hdf5_writer<double > writer(file, "double_float");
    for(int i=0; i<5; ++i) {
        double x = 1.1*i;
        writer.append(x);
    }
    H5Fclose(file);
}

BOOST_AUTO_TEST_CASE(array1d)
{
    hid_t file = H5Fcreate("array1d.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const int dim = 6;
    MArray1d a(boost::extents[dim]);
    Hdf5_writer<MArray1d_ref > writer(file, "array1d");
    for(int i=0; i<5; ++i) {
        for(int j = 0; j< dim; ++j) {
            a[j] = 1.1*j + 10*i;
        }
        writer.append(a);
    }
    H5Fclose(file);
}

