#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/distributed_fft2d.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/foundation/math_constants.h"
#include <complex>
#include <cmath>

// set DBGPRINT to 1 to print values for tolerance failures
#define DBGPRINT 0
#if DBGPRINT
#include <iostream>
#endif

// define FAILME nonzero to force failures
#define FAILME 0

BOOST_GLOBAL_FIXTURE(MPI_fixture);

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
        shape(shape_struct.shape), commxx_sptr(new Commxx()),
        distributed_fft2d(shape, commxx_sptr)
    {
    }
    ~Fixture()
    {
    }
    Shape_struct shape_struct;
    std::vector<int > shape;
    Commxx_sptr commxx_sptr;
    Distributed_fft2d distributed_fft2d;
};

struct Shape_struct2
{
    Shape_struct2() :
        shape(2)
    {
      shape[0] = 16;
      shape[1] = 8;
    }
    ~Shape_struct2()
    {

    }
    std::vector<int > shape;
};

struct Fixture2
{
    Fixture2() :
        shape(shape_struct2.shape), commxx_sptr(new Commxx()),
        distributed_fft2d(shape, commxx_sptr)
    {
    }
    ~Fixture2()
    {
    }
    Shape_struct2 shape_struct2;
    std::vector<int > shape;
    Commxx_sptr commxx_sptr;
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

BOOST_FIXTURE_TEST_CASE(transform_realtest, Fixture2)
{
    int lower = distributed_fft2d.get_lower();
    int upper = distributed_fft2d.get_upper();

    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc orig(
            boost::extents[extent_range(lower, upper)][shape[1]]);
    MArray2dc dest(
            boost::extents[extent_range(lower, upper)][shape[1]]);

    // the location of the constructed frequency space peak
    const int loc0 = 7;
    const int loc1 = 2;

    // The frequency corresponding to the array position
    const double freq0 = double(loc0)/shape[0];
    const double freq1 = double(loc1)/shape[1];

    const double dx0 = 1.0;
    const double dx1 = 1.0;

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
	// will with exp(2*pi*i*freq0*x0)*exp(2*pi*i*freq1*x1) which is
	// a delta function in frequency space
	orig[i0][i1] =
	  std::complex<double >( (cos(2.0*mconstants::pi*freq0*x0)*cos(2.0*mconstants::pi*freq1*x1) -
			 sin(2.0*mconstants::pi*freq0*x0)*sin(2.0*mconstants::pi*freq1*x1)) ,
			(sin(2.0*mconstants::pi*freq0*x0)*cos(2.0*mconstants::pi*freq1*x1) +
			 cos(2.0*mconstants::pi*freq0*x0)*sin(2.0*mconstants::pi*freq1*x1)) );
#if FAILME
	orig[i0][i1] *= 1.1;
#endif
	x1 += dx1;
      }
      x0 += dx0;
    }

    distributed_fft2d.transform(orig, dest);

    double norm = double(shape[0] * shape[1]);

    // I see deviations of a few *1e12
    double fft_tolerance = 1.0e-11;

    // All the numbers in the result should have a negligible complex part
    for (int i0=lower; i0<upper; ++i0) {
      for (int i1=0; i1<shape[1]; ++i1) {
	  BOOST_CHECK(std::abs(dest[i0][i1].imag()) < fft_tolerance);
#if DBGPRINT
	  if (std::abs(dest[i0][i1].imag()) >= fft_tolerance) {
	    std::cout << "imaginary part non-zero failure at (" <<
	      i0 << "," << i1 << "): " << std::abs(dest[i0][i1].imag()) << std::endl;
	  }
#endif // DBGPRINT
      }
    }

    // check the real parts.  All of them should be negligible except for
    // the one corresponding to the delta function location.
    for (int i0=lower; i0<upper; ++i0) {
      for (int i1=0; i1<shape[1]; ++i1) {
	if ( (i0==loc0)&&(i1==loc1) ) {
	  BOOST_CHECK(std::abs(dest[i0][i1].real()-norm) < fft_tolerance);
#if DBGPRINT
	  if (std::abs(dest[i0][i1].real()-norm) >= fft_tolerance) {
	    std::cout << "real part not correct failure at (" <<
	      i0 << "," << i1 << "): " << dest[i0][i1].real() << std::endl;
	  }
#endif // DBGPRINT
	} else {
	  BOOST_CHECK(std::abs(dest[i0][i1].real()) < fft_tolerance);
#if DBGPRINT
	  if (std::abs(dest[i0][i1].real()) >= fft_tolerance) {
	    std::cout << "real part not negligible failure at (" <<
	      i0 << "," << i1 << "): " << dest[i0][i1].real() << std::endl;
	  }
#endif // DBGPRINT
	}
      }
    }
}
