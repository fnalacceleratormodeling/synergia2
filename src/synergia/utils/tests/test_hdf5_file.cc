#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"

BOOST_AUTO_TEST_CASE(construct)
{
    Hdf5_file file("hdf5_file_empty.h5");
}

BOOST_AUTO_TEST_CASE(write_data)
{
    Hdf5_file file("hdf5_file_write.h5");
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
