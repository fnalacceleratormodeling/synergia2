#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/multi_array_offsets.h"

double
fn(int i, int j, int k)
{
    return 100 * i + 10 * j + k;
}

BOOST_AUTO_TEST_CASE(simple)
{
    MArray3d array(boost::extents[3][4][5]);
    for (int i = 0; i < array.shape()[0]; ++i) {
        for (int j = 0; j < array.shape()[1]; ++j) {
            for (int k = 0; k < array.shape()[2]; ++k) {
                array[i][j][k] = fn(i, j, k);
            }
        }
    }
    double * pointer;
    int i, j, k;

    i = 2;
    j = 0;
    k = 0;
    pointer = multi_array_offset(array, i, j, k);
    BOOST_CHECK_EQUAL(*pointer, fn(i,j,k));

    i = 0;
    j = 2;
    k = 0;
    pointer = multi_array_offset(array, i, j, k);
    BOOST_CHECK_EQUAL(*pointer, fn(i,j,k));

    i = 0;
    j = 0;
    k = 2;
    pointer = multi_array_offset(array, i, j, k);
    BOOST_CHECK_EQUAL(*pointer, fn(i,j,k));

    i = 2;
    j = 3;
    k = 4;
    pointer = multi_array_offset(array, i, j, k);
    BOOST_CHECK_EQUAL(*pointer, fn(i,j,k));
}

double
fn_2d(int i, int j)
{
    return 10 * i + j;
}

BOOST_AUTO_TEST_CASE(simple_2d)
{   
    MArray2d array(boost::extents[3][4]); 
    for (int i = 0; i < array.shape()[0]; ++i) {
        for (int j = 0; j < array.shape()[1]; ++j) {
                array[i][j] = fn_2d(i, j);
        }
    }
    double * pointer;
    int i, j, k;

    i = 2;
    j = 0;
    pointer = multi_array_offset(array, i, j);
    BOOST_CHECK_EQUAL(*pointer, fn_2d(i,j));

    i = 0;
    j = 2;
    pointer = multi_array_offset(array, i, j);
    BOOST_CHECK_EQUAL(*pointer, fn_2d(i,j));

    i = 2;
    j = 3;
    pointer = multi_array_offset(array, i, j);
    BOOST_CHECK_EQUAL(*pointer, fn_2d(i,j));
}
