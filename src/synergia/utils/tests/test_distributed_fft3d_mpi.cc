#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/distributed_fft3d.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/multi_array_print.h"
#include <complex>
#include <cmath>

// set DBGPRINT to 1 to print values for tolerance failures
#define DBGPRINT 0
#if DBGPRINT
#include <iostream>
#endif

// define FAILME to force failures
#define FAILME 0

BOOST_GLOBAL_FIXTURE(MPI_fixture)

// n.b. We use 0,1,2 here instead of x,y,z because
// we may use z,y,x ordering of arrays.
const int shape0 = 16;
const int shape1 = 4;
const int shape2 = 4;

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

struct Shape_struct2
{
    Shape_struct2() :
        shape(3)
    {
        shape[0] = 32;
        shape[1] = 16;
        shape[2] = 8;
    }
    ~Shape_struct2()
    {

    }
    std::vector<int > shape;
};

struct Fixture2
{
    Fixture2() :
        shape(shape_struct2.shape), distributed_fft3d(shape, Commxx(
                MPI_COMM_WORLD))
    {
    }
    ~Fixture2()
    {
    }
    Shape_struct2 shape_struct2;
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

BOOST_FIXTURE_TEST_CASE(transform_realtest, Fixture2)
{
    int lower = distributed_fft3d.get_lower();
    int upper = distributed_fft3d.get_upper();

    // rshape is the shape of the complex array but in numbers of real numbers
    // the array must be that large to accomodate the in-place transform
    std::vector<int > rshape(distributed_fft3d.get_padded_shape_real());
    MArray3d orig(
            boost::extents[extent_range(lower, upper)][rshape[1]][rshape[2]]);
    // cshape is the shape of the array in complex numbers
    std::vector<int > cshape(distributed_fft3d.get_padded_shape_complex());
    MArray3dc carray(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);

    // the location of the constructed frequency space peak
    const int loc0 = 13;
    const int loc1 = 7;
    const int loc2 = 2;

    // The frequency corresponding to the array position
    const double freq0 = double(loc0)/shape[0];
    const double freq1 = double(loc1)/shape[1];
    const double freq2 = double(loc2)/shape[2];
    
    const double dx0 = 1.0;
    const double dx1 = 1.0;
    const double dx2 = 1.0;

    // fill the input array
    double x0 = double(lower)*dx0;
    for (int i0=lower; i0<upper; ++i0) {
      if (i0 >= shape[0]/2) {
	x0 -= shape[0]*dx0;
      }
      double x1 = 0.0;
      for (int i1=0; i1<shape[1]; ++i1) {
	if (i1 >= shape[1]/2) {
	  x1 -= shape[1]*dx1;
	}
	double x2 = 0.0;
	for (int i2=0; i2<shape[2]; ++i2) {
	  if (i2 >= shape[2]/2) {
	    x2 -= shape[2]*dx2;
	  }
	  orig[i0][i1][i2] =
	    cos(2.0*M_PI*freq0*x0) *
	    cos(2.0*M_PI*freq1*x1) *
	    cos(2.0*M_PI*freq2*x2);
#if FAILME
	  orig[i0][i1][i2] *= 1.1;
	  if (i0 > shape[0]/2) {
	    orig[i0][i1][i2] = 0.0;
	  }
#endif
	  x2 += dx2;
	}
	x1 += dx1;
      }
      x0 += dx0;
    }

    distributed_fft3d.transform(orig, carray);
    
    // it's over 8 because the single frequency spike is reflected 8 times, even
    // though we only see four of them (because we're looking at the
    // half-complex transform result.)
    double norm = double(shape[0] * shape[1] * shape[2]/8);
    
    // I see deviations of a few *1e12
    double fft_tolerance = 1.0e-11;
    
    // All the numbers in the result should have a negligible complex part
    for (int i0=lower; i0<upper; ++i0) {
      for (int i1=0; i1<rshape[1]; ++i1) {
	for (int i2=0; i2<rshape[2]; ++i2) {
	  BOOST_CHECK(std::abs(carray[i0][i1][i2].imag()) < fft_tolerance);
#if DBGPRINT
	  if (std::abs(carray[i0][i1][i2].imag()) >= fft_tolerance) {
	    std::cout << "imaginary part non-zero failure at (" <<
	      i0 << "," << i1 << "," << i2 << "): " << std::abs(carray[i0][i1][i2].imag()) << std::endl;
	  }
#endif // DBGPRINT
	}
      }
    }

    // check the real parts.  All of them should be negligible except for
    // the ones corresponding to the location and the mirrors.
    for (int i0=lower; i0<upper; ++i0) {
      for (int i1=0; i1<cshape[1]; ++i1) {
	for (int i2=0; i2<cshape[2]; ++i2) {
	  if ( (i0==loc0)&&(i1==loc1)&&(i2==loc2) ||
	       (i0==loc0)&&(i1==(shape[1]-loc1))&&(i2==loc2) ||
	       (i0==(shape[0]-loc0))&&(i1==loc1)&&(i2==loc2) ||
	       (i0==(shape[0]-loc0))&&(i1==(shape[1]-loc1))&&(i2==loc2) ) {
	    BOOST_CHECK(std::abs(carray[i0][i1][i2].real()-norm) < fft_tolerance);
#if DBGPRINT
	    if (std::abs(carray[i0][i1][i2].real()-norm) >= fft_tolerance) {
	      std::cout << "real part not correct failure at (" <<
		i0 << "," << i1 << "," << i2 << "): " << carray[i0][i1][i2].real() << std::endl;
	    }
#endif // DBGPRINT
	  } else {
	    BOOST_CHECK(std::abs(carray[i0][i1][i2].real()) < fft_tolerance);
#if DBGPRINT
	    if (std::abs(carray[i0][i1][i2].real()) >= fft_tolerance) {
	      std::cout << "real part not negligible failure at (" <<
		i0 << "," << i1 << "," << i2 << "): " << carray[i0][i1][i2].real() << std::endl;
	    }
#endif // DBGPRINT
	  }
	}
      }
    }
}

