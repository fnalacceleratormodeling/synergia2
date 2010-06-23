#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "utils/distributed_fft3d.h"
#include "utils/boost_test_mpi_fixture.h"
#include "utils/multi_array_check_equal.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

// n.b. We use 0,1,2 here instead of x,y,z because
// we may use z,y,x ordering of arrays.
const int shape0 = 3;
const int shape1 = 4;
const int shape2 = 5;

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

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
}

BOOST_FIXTURE_TEST_CASE(construct_measure, Fixture)
{
    Distributed_fft3d distributed_fft3d_measure(shape, Commxx(MPI_COMM_WORLD),
            FFTW_MEASURE);
}

BOOST_FIXTURE_TEST_CASE(get_lower, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft3d.get_lower(), 0);
}

BOOST_FIXTURE_TEST_CASE(get_upper, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft3d.get_upper(), shape0);
}

BOOST_FIXTURE_TEST_CASE(get_guard_lower, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft3d.get_guard_lower(), 0);
}

BOOST_FIXTURE_TEST_CASE(get_guard_upper, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft3d.get_guard_upper(), shape0);
}

BOOST_FIXTURE_TEST_CASE(get_offset, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft3d.get_offset(), 0);
}

BOOST_FIXTURE_TEST_CASE(get_local_size, Fixture)
{
    BOOST_CHECK(distributed_fft3d.get_local_size() >= shape0*shape1*shape2);
}

BOOST_FIXTURE_TEST_CASE(get_padded_shape_real, Fixture)
{
    std::vector<int > got_shape(distributed_fft3d.get_padded_shape_real());
    BOOST_CHECK_EQUAL(got_shape[0], shape[0]);
    BOOST_CHECK_EQUAL(got_shape[1], shape[1]);
    BOOST_CHECK_EQUAL(got_shape[2], 2 * (shape[2] / 2 + 1));
}

BOOST_FIXTURE_TEST_CASE(get_padded_shape_complex, Fixture)
{
    std::vector<int > got_shape(distributed_fft3d.get_padded_shape_complex());
    BOOST_CHECK_EQUAL(got_shape[0], shape[0]);
    BOOST_CHECK_EQUAL(got_shape[1], shape[1]);
    BOOST_CHECK_EQUAL(got_shape[2], shape[2] / 2 + 1);
}

const double tolerance = 1.0e-12;
BOOST_FIXTURE_TEST_CASE(transform_roundtrip, Fixture)
{
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d rarray(boost::extents[rshape[0]][rshape[1]][rshape[2]]);
    MArray3d orig(boost::extents[rshape[0]][rshape[1]][rshape[2]]);
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(boost::extents[cshape[0]][cshape[1]][cshape[2]]);
    for (int j = 0; j < shape[0]; ++j) {
        for (int k = 0; k < shape[1]; ++k) {
            for (int l = 0; l < shape[2]; ++l) {
                rarray[j][k][l] = 1.1 * l + 10 * k + 100 * j;
                orig[j][k][l] = 1.1 * l + 10 * k + 100 * j;
            }
        }
    }

    distributed_fft3d.transform(rarray, carray);
    distributed_fft3d.inv_transform(carray, rarray);
    double norm = shape[0] * shape[1] * shape[2];
    for (int j = 0; j < shape[0]; ++j) {
        for (int k = 0; k < shape[1]; ++k) {
            for (int l = 0; l < shape[2]; ++l) {
                rarray[j][k][l] *= 1.0 / norm;
            }
        }
    }

    // zero out padded region
    for (int j = 0; j < rshape[0]; ++j) {
        for (int k = 0; k < rshape[1]; ++k) {
            for (int l = shape[2]; l < rshape[2]; ++l) {
                rarray[j][k][l] = 0.0;
                orig[j][k][l] = 0.0;
            }
        }
    }
    multi_array_check_equal(orig, rarray, tolerance);
}
