#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/distributed_fft3d.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/multi_array_print.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

// n.b. We use 0,1,2 here instead of x,y,z because
// we may use z,y,x ordering of arrays.
const int shape0 = 4;
const int shape1 = 2;
const int shape2 = 2;

struct Shape_struct
{
    Shape_struct() :
        shape(3)
    {
        shape[0] = shape0;
        shape[1] = shape1;
        shape[2] = shape2;
    }
    ~Shape_struct()
    {

    }
    std::vector<int > shape;
};

struct Fixture
{
    Fixture() :
        shape(shape_struct.shape), distributed_fft3d(shape, Commxx(
                MPI_COMM_WORLD))
    {
    }
    ~Fixture()
    {
    }
    Shape_struct shape_struct;
    std::vector<int > shape;
    Distributed_fft3d distributed_fft3d;
};

BOOST_FIXTURE_TEST_CASE(get_uppers, Fixture)
{
    std::vector<int > got_uppers(distributed_fft3d.get_uppers());
    int size = distributed_fft3d.get_comm().get_size();
    BOOST_CHECK_EQUAL(got_uppers.size(), size);
    BOOST_CHECK_EQUAL(got_uppers[size-1], shape[0]);
}

BOOST_FIXTURE_TEST_CASE(get_lengths, Fixture)
{
    std::vector<int > got_lengths(distributed_fft3d.get_lengths());
    int size = distributed_fft3d.get_comm().get_size();
    BOOST_CHECK_EQUAL(got_lengths.size(), size);
    int total_length = 0;
    for (int i = 0; i < size; ++i) {
        total_length += got_lengths[i];
    }
    BOOST_CHECK_EQUAL(total_length, shape[0]*shape[1]*shape[2]);
}

const double tolerance = 1.0e-10;
BOOST_FIXTURE_TEST_CASE(transform_roundtrip, Fixture)
{
    int lower = distributed_fft3d.get_lower();
    int upper = distributed_fft3d.get_upper();
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d rarray(
            boost::extents[extent_range(lower, upper)][rshape[1]][rshape[2]]);
    MArray3d orig(
            boost::extents[extent_range(lower, upper)][rshape[1]][rshape[2]]);
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
            for (int k = 0; k < shape[2]; ++k) {
                rarray[i][j][k] = 1.1 * k + 10 * j + 100 * i;
                orig[i][j][k] = 1.1 * k + 10 * j + 100 * i;
            }
        }
    }

    distributed_fft3d.transform(rarray, carray);
    distributed_fft3d.inv_transform(carray, rarray);
    double norm = shape[0] * shape[1] * shape[2];
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
            for (int k = 0; k < shape[2]; ++k) {
                rarray[i][j][k] *= 1.0 / norm;
            }
        }
    }

    // zero out padded region
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < rshape[1]; ++j) {
            for (int k = shape[2]; k < rshape[2]; ++k) {
                rarray[i][j][k] = 0.0;
                orig[i][j][k] = 0.0;
            }
        }
    }
    multi_array_check_equal(orig, rarray, tolerance);
}

BOOST_FIXTURE_TEST_CASE(transform_bad_in_offset, Fixture)
{
    int lower = distributed_fft3d.get_lower();
    int upper = distributed_fft3d.get_upper();
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d rarray(
            boost::extents[extent_range(lower-1, upper-1)][rshape[1]][rshape[2]]);
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);

    bool caught_error = false;
    try {
        distributed_fft3d.transform(rarray, carray);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(transform_bad_out_offset, Fixture)
{
    int lower = distributed_fft3d.get_lower();
    int upper = distributed_fft3d.get_upper();
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d rarray(
            boost::extents[extent_range(lower, upper)][rshape[1]][rshape[2]]);
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(
            boost::extents[extent_range(lower-1, upper-1)][cshape[1]][cshape[2]]);

    bool caught_error = false;
    try {
        distributed_fft3d.transform(rarray, carray);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_bad_in_offset, Fixture)
{
    int lower = distributed_fft3d.get_lower();
    int upper = distributed_fft3d.get_upper();
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d rarray(
            boost::extents[extent_range(lower, upper)][rshape[1]][rshape[2]]);
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(
            boost::extents[extent_range(lower-1, upper-1)][cshape[1]][cshape[2]]);

    bool caught_error = false;
    try {
        distributed_fft3d.inv_transform(carray, rarray);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_bad_out_offset, Fixture)
{
    int lower = distributed_fft3d.get_lower();
    int upper = distributed_fft3d.get_upper();
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d rarray(
            boost::extents[extent_range(lower-1, upper-1)][rshape[1]][rshape[2]]);
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);

    bool caught_error = false;
    try {
        distributed_fft3d.inv_transform(carray, rarray);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

