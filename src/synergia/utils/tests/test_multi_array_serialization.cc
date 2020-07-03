#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_serialization.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/serialization_files.h"

#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-14;

BOOST_AUTO_TEST_CASE(test1d)
{
    MArray1d a(boost::extents[7]);
    for (unsigned int i = 0; i < a.shape()[0]; ++i) {
        a[i] = 1.1 * i;
    }
    xml_save<MArray1d > (a, "multi_array1d.xml");
    MArray1d b;
    xml_load<MArray1d > (b, "multi_array1d.xml");
    multi_array_check_equal(a, b, tolerance);
}

BOOST_AUTO_TEST_CASE(test2d)
{
    MArray2d a(boost::extents[7][9]);
    for (unsigned int i = 0; i < a.shape()[0]; ++i) {
        for (unsigned int j = 0; j < a.shape()[1]; ++j) {
            a[i][j] = 11.1 * i + 0.93 * j;
        }
    }
    xml_save<MArray2d > (a, "multi_array2d.xml");
    MArray2d b;
    xml_load<MArray2d > (b, "multi_array2d.xml");
    multi_array_check_equal(a, b, tolerance);
}

BOOST_AUTO_TEST_CASE(test2d_b)
{
    MArray2d a(boost::extents[6][6]);
    for (unsigned int i = 0; i < a.shape()[0]; ++i) {
        for (unsigned int j = 0; j < a.shape()[1]; ++j) {
            a[i][j] = 0.0;
        }
    }
    for (unsigned int i=0; i < a.shape()[0]; ++i) {
        a[i][i] = 1.0*i;
        if (i%2 == 0) {
            a[i][i+1] = i*1.0+0.5;
        } else {
            a[i][i-1] = i*1.0-0.5;
        }
    }
    std::cout << "saving array:" << std::endl;
    for (int i=0; i<6; ++i) {
        for (int j=0; j<6; ++j) {
            if (j != 0) {
                std::cout << " ";
            }
            std::cout << a[i][j];
        }
        std::cout << std::endl;
    }
    xml_save(a, "multi_array2d_b.xml");
    MArray2d b;
    xml_load(b, "multi_array2d_b.xml");
    std::cout << "read back" << std::endl;
    for (int i=0; i<6; ++i) {
        for (int j=0; j<6; ++j) {
            if (j != 0) {
                std::cout << " ";
            }
            std::cout << b[i][j];
        }
        std::cout << std::endl;
    }
    multi_array_check_equal(a, b, tolerance);
}
