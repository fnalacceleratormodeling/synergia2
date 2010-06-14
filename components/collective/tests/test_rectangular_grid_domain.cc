#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/collective/rectangular_grid_domain.h"
#include "utils/boost_test_mpi_fixture.h"
//BOOST_GLOBAL_FIXTURE(MPI_fixture)
//;

const double tolerance = 1.0e-12;

const double domain_min = -1.0;
const double domain_max = 1.0;
const double domain_offset = 5.0;
const int grid_size = 2;
const bool is_periodic = false;

struct Rgd_fixture
{
    Rgd_fixture() :
        physical_size(3), physical_offset(3), grid_shape(3)
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
            grid_shape[i] = grid_size;
        }
        rectangular_grid_domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(physical_size, physical_offset,
                        grid_shape, is_periodic));
    }

    ~Rgd_fixture()
    {
    }

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Rectangular_grid_domain_sptr rectangular_grid_domain_sptr;
};

BOOST_FIXTURE_TEST_CASE(construct, Rgd_fixture)
{
}

BOOST_FIXTURE_TEST_CASE(get_physical_size, Rgd_fixture)
{
    std::vector<double > physical_size =
            rectangular_grid_domain_sptr->get_physical_size();
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(physical_size.at(i), domain_max - domain_min, tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_physical_offset, Rgd_fixture)
{
    std::vector<double > physical_offset =
            rectangular_grid_domain_sptr->get_physical_offset();
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(physical_offset.at(i), domain_offset, tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(get_grid_shape, Rgd_fixture)
{
    std::vector<int > grid_shape =
            rectangular_grid_domain_sptr->get_grid_shape();
    for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_CLOSE(grid_shape.at(i), grid_size, tolerance);
    }
}

BOOST_FIXTURE_TEST_CASE(is_periodic_, Rgd_fixture)
{
    BOOST_CHECK_EQUAL(rectangular_grid_domain_sptr->is_periodic(), is_periodic);
}

BOOST_FIXTURE_TEST_CASE(get_leftmost_indices_offsets, Rgd_fixture)
{
    int ix, iy, iz;
    double offx, offy, offz;
    double testx = domain_offset - 0.1, testy = domain_offset + 0.2, testz =
            domain_offset + 0.9;
    rectangular_grid_domain_sptr->get_leftmost_indices_offsets(testx, testy,
            testz, ix, iy, iz, offx, offy, offz);
    BOOST_CHECK_EQUAL(ix,0);
    BOOST_CHECK_CLOSE(offx, 0.9, tolerance);
    BOOST_CHECK_EQUAL(iy,1);
    BOOST_CHECK_CLOSE(offy, 0.2, tolerance);
    BOOST_CHECK_EQUAL(iz,1);
    BOOST_CHECK_CLOSE(offz, 0.9, tolerance);
}
