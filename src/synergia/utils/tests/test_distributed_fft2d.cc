#include "synergia/utils/catch.hpp"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/distributed_fft2d.h"

#include <complex>
#include <cmath>

// set DBGPRINT to 1 to print values for tolerance failures
#define DBGPRINT 0

#if DBGPRINT
#include <iostream>
#endif

// define FAILME nonzero to force failures
#define FAILME 0

const int shape0 = 16;
const int shape1 = 8;

TEST_CASE("construct")
{
    REQUIRE_NOTHROW(Distributed_fft2d());
}


TEST_CASE("methods")
{
    // construct the fft2d object
    Distributed_fft2d fft;

    // construct the communicator so that every rank does
    // the full calculation, by dividing the world into 
    // subgroups of size 1
    auto comm_world = std::make_shared<Commxx>(Commxx::World);
    auto comm = comm_world->divide(1);

    fft.construct({shape0, shape1}, comm);

    SECTION("roundtrip normalization")
    {
        double normalization = fft.get_roundtrip_normalization();
        CHECK(normalization == 1.0/(shape0*shape1));
    }

    SECTION("get shape")
    {
        auto s = fft.get_shape();

        CHECK(s.size() == 2);
        CHECK(s[0] == shape0);
        CHECK(s[1] == shape1);
    }

    SECTION("get comm")
    {
        CHECK(fft.get_comm() == comm);
    }

    SECTION("get lower")
    {
        CHECK(fft.get_lower() == 0);
        CHECK(fft.get_upper() == shape0);
    }
}


const double tolerance = 1.0e-12;
TEST_CASE("transform_roundtrip")
{
    // construct the fft2d object
    Distributed_fft2d fft;

    // construct the communicator so that every rank does
    // the full calculation, by dividing the world into 
    // subgroups of size 1
    auto comm_world = std::make_shared<Commxx>(Commxx::World);

#ifdef Kokkos_ENABLE_CUDA
    int subgroup_size = 1;
#else
    int subgroup_size = Commxx::world_size();
#endif

    auto comm = comm_world->divide(subgroup_size);

    fft.construct({shape0, shape1}, comm);

    karray1d src("src", shape0*shape1*2);
    karray1d orig("orig", shape0*shape1*2);

    for(int j=0; j<shape0; ++j)
    {
        for(int k=0; k<shape1; ++k)
        {
            src(j*shape1*2+k*2+0) = 1.1*k+10*j;
            src(j*shape1*2+k*2+1) = 1.2*k+10*j;

            orig(j*shape1*2+k*2+0) = 1.1*k+10*j;
            orig(j*shape1*2+k*2+1) = 1.2*k+10*j;
        }
    }

    karray1d_dev d_src = Kokkos::create_mirror_view(src);
    Kokkos::deep_copy(d_src, src);

    karray1d_dev d_dst("dest", shape0*shape1*2);

    fft.transform(d_src, d_dst);
    fft.inv_transform(d_dst, d_src);

    Kokkos::deep_copy(src, d_src);

    auto norm = fft.get_roundtrip_normalization();
    int lower = fft.get_lower();
    int upper = fft.get_upper();

    for(int j=lower; j<upper; ++j)
    {
        for(int k=0; k<shape1; ++k)
        {
            int idx_real = j*shape1*2+k*2+0;
            int idx_imag = j*shape1*2+k*2+1;

            CHECK(src(idx_real)*norm == Approx(orig(idx_real)).margin(tolerance));
            CHECK(src(idx_imag)*norm == Approx(orig(idx_imag)).margin(tolerance));
        }
    }
}

TEST_CASE("transform_realtest")
{
    // construct the fft2d object
    Distributed_fft2d fft;

    // construct the communicator so that every rank does
    // the full calculation, by dividing the world into 
    // subgroups of size 1
    auto comm_world = std::make_shared<Commxx>(Commxx::World);

#ifdef Kokkos_ENABLE_CUDA
    int subgroup_size = 1;
#else
    int subgroup_size = Commxx::world_size();
#endif

    auto comm = comm_world->divide(subgroup_size);

    fft.construct({shape0, shape1}, comm);

    karray1d orig("orig", shape0*shape1*2);
    karray1d dest("dest", shape0*shape1*2);

    karray1d_dev d_orig = Kokkos::create_mirror_view(orig);
    karray1d_dev d_dest = Kokkos::create_mirror_view(dest);

    int lower = fft.get_lower();
    int upper = fft.get_upper();

    // the location of the constructed frequency space peak
    const int loc0 = 7;
    const int loc1 = 2;

    // in range
    REQUIRE(loc0 > 0);
    REQUIRE(loc0 < shape0);

    REQUIRE(loc1 > 0);
    REQUIRE(loc1 < shape1);

    // The frequency corresponding to the array position
    const double freq0 = double(loc0)/shape0;
    const double freq1 = double(loc1)/shape1;

    const double dx0 = 1.0;
    const double dx1 = 1.0;

    // fill the input array
    double x0 = lower * dx0;

    const double pi = mconstants::pi;
    
    for (int i0=lower; i0<upper; ++i0) 
    {
        if (i0 >= shape0/2) x0 -= shape0*dx0;

        double x1 = 0.0; 

        for (int i1=0; i1<shape1; ++i1) 
        { 
            if (i1 >= shape1/2) x1 -= shape1*dx1;

            // will with exp(2*pi*i*freq0*x0)*exp(2*pi*i*freq1*x1) which is
            // a delta function in frequency space
            double real = (cos(2.0*pi*freq0*x0)*cos(2.0*pi*freq1*x1) 
                    - sin(2.0*pi*freq0*x0)*sin(2.0*pi*freq1*x1));

            double imag = (sin(2.0*pi*freq0*x0)*cos(2.0*pi*freq1*x1) 
                    + cos(2.0*pi*freq0*x0)*sin(2.0*pi*freq1*x1));

            orig(i0*shape1*2+i1*2+0) = real;
            orig(i0*shape1*2+i1*2+1) = imag;

#if FAILME
            orig(i0*shape1*2+i1*2+0) *= 1.1;
            orig(i0*shape1*2+i1*2+1) *= 1.1;
#endif

            x1 += dx1; 
        } 

        x0 += dx0;
    }

    // copy orig to device
    Kokkos::deep_copy(d_orig, orig);

    // transform
    fft.transform(d_orig, d_dest);

    // and copy the dest back to host
    Kokkos::deep_copy(dest, d_dest);

    // now verify the results
    double norm = double(shape0 * shape1);

    // I see deviations of a few *1e12
    double fft_tolerance = 1.0e-11;

    // All the numbers in the result should have a negligible complex part
    for (int i0=lower; i0<upper; ++i0) 
    { 
        for (int i1=0; i1<shape1; ++i1) 
        { 
            double imag = dest(i0*shape1*2 + i1*2 + 1);
            CHECK(std::abs(imag) < fft_tolerance);

#if DBGPRINT
            if (std::abs(imag) >= fft_tolerance) 
            {
                std::cout 
                    << "imaginary part non-zero failure at (" 
                    << i0 << "," << i1 << "): " 
                    << imag << std::endl; 
            }
#endif // DBGPRINT 
        } 
    }

    // check the real parts.  All of them should be negligible except for
    // the one corresponding to the delta function location.
    for (int i0=lower; i0<upper; ++i0) 
    { 
        for (int i1=0; i1<shape1; ++i1) 
        { 
            double real = dest(i0*shape1*2 + i1*2 + 0);

            if ( i0==loc0 && i1==loc1 ) 
            { 
                CHECK(std::abs(real-norm) < fft_tolerance);

#if DBGPRINT
                if (std::abs(real-norm) >= fft_tolerance) 
                { 
                    std::cout 
                        << "real part not correct failure at (" 
                        << i0 << "," << i1 << "): " 
                        << real << std::endl; 
                }
#endif // DBGPRINT 
            } 
            else 
            { 
                CHECK(std::abs(real) < fft_tolerance);

#if DBGPRINT
                if (std::abs(real) >= fft_tolerance) 
                { 
                    std::cout 
                        << "real part not negligible failure at (" 
                        << i0 << "," << i1 << "): " 
                        << real << std::endl; 
                }
#endif // DBGPRINT 
            } 
        } 
    }
}


#if 0

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
}

BOOST_FIXTURE_TEST_CASE(get_comm, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft2d.get_comm().get(), MPI_COMM_WORLD);
}

BOOST_FIXTURE_TEST_CASE(get_lower, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft2d.get_lower(), 0);
}

BOOST_FIXTURE_TEST_CASE(get_upper, Fixture)
{
    BOOST_CHECK_EQUAL(distributed_fft2d.get_upper(), shape0);
}

BOOST_FIXTURE_TEST_CASE(get_uppers, Fixture)
{
    std::vector<int > got_uppers(distributed_fft2d.get_uppers());
    BOOST_CHECK_EQUAL(got_uppers.size(), 1);
    BOOST_CHECK_EQUAL(got_uppers[0], shape[0]);
}

BOOST_FIXTURE_TEST_CASE(get_lengths, Fixture)
{
    std::vector<int > got_lengths(distributed_fft2d.get_lengths());
    BOOST_CHECK_EQUAL(got_lengths.size(), 1);
    BOOST_CHECK_EQUAL(got_lengths[0], shape[0]*shape[1]);
}

BOOST_FIXTURE_TEST_CASE(get_shape, Fixture)
{
    std::vector<int > got_shape(distributed_fft2d.get_shape());
    BOOST_CHECK_EQUAL(got_shape[0], shape[0]);
    BOOST_CHECK_EQUAL(got_shape[1], shape[1]);
}

BOOST_FIXTURE_TEST_CASE(get_roundtrip_normalization, Fixture)
{
    double normalization = distributed_fft2d.get_roundtrip_normalization();
    BOOST_CHECK_EQUAL(normalization, 1.0/(shape[0]*shape[1]));
}

BOOST_FIXTURE_TEST_CASE(data_size,Fixture)
{
    int data_in = distributed_fft2d.get_data_size();
    int data_out = distributed_fft2d.get_workspace_size();
    //std::cout << "data_in = " << data_in << std::endl;
    //std::cout << "data_out = " << data_out << std::endl;
    BOOST_CHECK_EQUAL(data_in, data_out);
}

const double tolerance = 1.0e-12;
BOOST_FIXTURE_TEST_CASE(transform_roundtrip, Fixture)
{
    std::vector<int > grid_shape(distributed_fft2d.get_shape());
    MArray2dc src_array(boost::extents[grid_shape[0]][grid_shape[1]]);
    MArray2dc orig_array(boost::extents[grid_shape[0]][grid_shape[1]]);
    MArray2dc dest_array(boost::extents[grid_shape[0]][grid_shape[1]]);
    for (int j = 0; j < shape[0]; ++j) {
        for (int k = 0; k < shape[1]; ++k) {
            src_array[j][k] = std::complex<double > (1.1 * k + 10 * j, 1.2 * k + 10 * j);
            orig_array[j][k] = std::complex<double > (1.1 * k + 10 * j, 1.2 * k + 10 * j);
        }
    }

    distributed_fft2d.transform(src_array, dest_array);
    distributed_fft2d.inv_transform(dest_array, src_array);
    for (int j = 0; j < shape[0]; ++j) {
        for (int k = 0; k < shape[1]; ++k) {
            src_array[j][k]
                    *= distributed_fft2d.get_roundtrip_normalization();
        }
    }

    multi_complex_array_check_equal(orig_array, src_array, tolerance);
}

BOOST_FIXTURE_TEST_CASE(transform_padded_shape0, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    shape[0] += 1;
    MArray2dc src_array(
            boost::extents[extent_range(-1, shape[0])][shape[1]]);
    MArray2dc dest_array(boost::extents[shape[0]][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.transform(src_array, dest_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == false);
}

BOOST_FIXTURE_TEST_CASE(transform_bad_in_shape, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    shape[1] += 1;
    MArray2dc src_array(boost::extents[shape[0]][shape[1]]);
    MArray2dc dest_array(boost::extents[shape[0]][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.transform(src_array, dest_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(transform_padded_out_shape0, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(boost::extents[shape[0]][shape[1]]);
    shape[0] += 1;
    MArray2dc dest_array(
            boost::extents[extent_range(-1, shape[0])][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.transform(src_array, dest_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == false);
}

BOOST_FIXTURE_TEST_CASE(transform_bad_out_shape, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(boost::extents[shape[0]][shape[1]]);
    shape[1] += 1;
    MArray2dc dest_array(boost::extents[shape[0]][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.transform(src_array, dest_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_padded_in_shape0, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(boost::extents[shape[0]][shape[1]]);
    shape[0] += 1;
    MArray2dc dest_array(
            boost::extents[extent_range(-1, shape[0])][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.inv_transform(dest_array, src_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == false);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_bad_in_shape, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    MArray2dc src_array(boost::extents[shape[0]][shape[1]]);
    shape[1] += 1;
    MArray2dc dest_array(boost::extents[shape[0]][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.inv_transform(dest_array, src_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == true);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_padded_out_shape0, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    shape[0] += 1;
    MArray2dc src_array(
            boost::extents[extent_range(-1, shape[0])][shape[1]]);
    MArray2dc dest_array(boost::extents[shape[0]][shape[1]]);

    bool caught_error = false;
    try {
        distributed_fft2d.inv_transform(dest_array, src_array);
    }
    catch (std::runtime_error) {
        caught_error = true;
    }

    BOOST_CHECK(caught_error == false);
}

BOOST_FIXTURE_TEST_CASE(inv_transform_bad_out_shape, Fixture)
{
    std::vector<int > shape(distributed_fft2d.get_shape());
    shape[1] += 1;
    MArray2dc src_array(boost::extents[shape[0]][shape[1]]);
    MArray2dc dest_array(boost::extents[shape[0]][shape[1]]);

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

#endif
