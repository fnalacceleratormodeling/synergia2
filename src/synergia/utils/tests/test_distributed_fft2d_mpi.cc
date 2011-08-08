#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/distributed_fft2d.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/multi_array_print.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

const int shape0 = 16;
const int shape1 = 4;

struct Shape_struct
{
    Shape_struct() :
        shape(2)
    {
        shape[0] = shape0;
        shape[1] = shape1;
    }
    ~Shape_struct()
    {

    }
    std::vector<int > shape;
};

struct Fixture
{
    Fixture() :
        shape(shape_struct.shape), distributed_fft2d(shape, Commxx(
                MPI_COMM_WORLD))
    {
    }
    ~Fixture()
    {
    }
    Shape_struct shape_struct;
    std::vector<int > shape;
    Distributed_fft2d distributed_fft2d;
};

BOOST_FIXTURE_TEST_CASE(get_uppers, Fixture)
{
    std::vector<int > got_uppers(distributed_fft2d.get_uppers());
    int size = distributed_fft2d.get_comm().get_size();
    BOOST_CHECK_EQUAL(got_uppers.size(), size);
    BOOST_CHECK_EQUAL(got_uppers[size-1], shape[0]);
}

BOOST_FIXTURE_TEST_CASE(get_lengths, Fixture)
{
    std::vector<int > got_lengths(distributed_fft2d.get_lengths());
    int size = distributed_fft2d.get_comm().get_size();
    BOOST_CHECK_EQUAL(got_lengths.size(), size);
    int total_length = 0;
    for (int i = 0; i < size; ++i) {
        total_length += got_lengths[i];
    }
    BOOST_CHECK_EQUAL(total_length, shape[0]*shape[1]);
}

const double tolerance = 1.0e-10;
BOOST_FIXTURE_TEST_CASE(transform_roundtrip, Fixture)
{
    int lower = distributed_fft2d.get_lower();
    int upper = distributed_fft2d.get_upper();
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);
    MArray2dc orig_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);
    MArray2dc dest_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
                src_array[i][j] = 1.1 * j + 10 * i;
                orig_array[i][j] = 1.1 * j + 10 * i;
        }
    }

    distributed_fft2d.transform(src_array, dest_array);
    distributed_fft2d.inv_transform(dest_array, src_array);
    double norm = shape[0] * shape[1];
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
                src_array[i][j] *= 1.0 / norm;
        }
    }

    multi_complex_array_check_equal(orig_array, src_array, tolerance);
}

BOOST_FIXTURE_TEST_CASE(transform_bad_in_offset, Fixture)
{
    int lower = distributed_fft2d.get_lower();
    int upper = distributed_fft2d.get_upper();
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(
            boost::extents[extent_range(lower-1, upper-1)][shape[1]]);
    MArray2dc dest_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.transform(src_array, dest_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(transform_bad_out_offset, Fixture)
{
    int lower = distributed_fft2d.get_lower();
    int upper = distributed_fft2d.get_upper();
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);
    MArray2dc dest_array(
            boost::extents[extent_range(lower-1, upper-1)][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.transform(src_array, dest_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_bad_in_offset, Fixture)
{
    int lower = distributed_fft2d.get_lower();
    int upper = distributed_fft2d.get_upper();
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);
    MArray2dc dest_array(
            boost::extents[extent_range(lower-1, upper-1)][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.inv_transform(dest_array, src_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_bad_out_offset, Fixture)
{
    int lower = distributed_fft2d.get_lower();
    int upper = distributed_fft2d.get_upper();
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(
            boost::extents[extent_range(lower-1, upper-1)][shape[1]]);
    MArray2dc dest_array(
            boost::extents[extent_range(lower, upper)][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.inv_transform(dest_array, src_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}
