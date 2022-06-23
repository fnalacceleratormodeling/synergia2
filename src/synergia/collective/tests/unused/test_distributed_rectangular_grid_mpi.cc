#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include "synergia/collective/distributed_rectangular_grid.h"
#include "rectangular_grid_domain_fixture.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
int grid_midpoint0 = grid_size0 / 2;

struct Distributed_rectangular_grid_domain_fixture
{
    Distributed_rectangular_grid_domain_fixture() :
        physical_size(3), physical_offset(3), grid_shape(3), comm_sptr(
                new Commxx), offsets(comm_sptr->get_size()), counts(
                comm_sptr->get_size())
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
        }
        grid_shape[0] = grid_size0;
        grid_shape[1] = grid_size1;
        grid_shape[2] = grid_size2;
        rectangular_grid_domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(physical_size, physical_offset,
                        grid_shape, false));
        decompose_1d(*comm_sptr, grid_shape[0], offsets, counts);
        int lower = offsets[comm_sptr->get_rank()];
        int upper = offsets[comm_sptr->get_rank()] + counts[comm_sptr->get_rank()];
        distributed_rectangular_grid_sptr = Distributed_rectangular_grid_sptr(
                new Distributed_rectangular_grid(rectangular_grid_domain_sptr,
                        lower, upper, comm_sptr));
    }

    ~Distributed_rectangular_grid_domain_fixture()
    {
    }

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Commxx_sptr comm_sptr;
    std::vector<int > offsets, counts;
    Rectangular_grid_domain_sptr rectangular_grid_domain_sptr;
    Distributed_rectangular_grid_sptr distributed_rectangular_grid_sptr;
};

struct Distributed_rectangular_grid_domain_with_empty_fixture
{
    Distributed_rectangular_grid_domain_with_empty_fixture() :
        physical_size(3), physical_offset(3), grid_shape(3), comm_sptr(
                new Commxx), offsets(comm_sptr->get_size()), counts(
                comm_sptr->get_size())
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
        }
        grid_shape[0] = grid_size0;
        grid_shape[1] = grid_size1;
        grid_shape[2] = grid_size2;
        rectangular_grid_domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(physical_size, physical_offset,
                        grid_shape, false));
        decompose_1d(*comm_sptr, 2 * grid_shape[0], offsets, counts);
        int lower = std::min(offsets[comm_sptr->get_rank()], grid_shape[0]);
        int upper = std::min(
                offsets[comm_sptr->get_rank()] + counts[comm_sptr->get_rank()],
                grid_shape[0]);
        distributed_rectangular_grid_sptr = Distributed_rectangular_grid_sptr(
                new Distributed_rectangular_grid(rectangular_grid_domain_sptr,
                        lower, upper, comm_sptr));
    }

    ~Distributed_rectangular_grid_domain_with_empty_fixture()
    {
    }

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Commxx_sptr comm_sptr;
    std::vector<int > offsets, counts;
    Rectangular_grid_domain_sptr rectangular_grid_domain_sptr;
    Distributed_rectangular_grid_sptr distributed_rectangular_grid_sptr;
};

// jfa: this is an ugly duplication in order to change one bool. I don't know how
// to do it more elegantly.
struct Distributed_rectangular_grid_domain_fixture_periodic
{
    Distributed_rectangular_grid_domain_fixture_periodic() :
        physical_size(3), physical_offset(3), grid_shape(3), comm_sptr(
                new Commxx), offsets(comm_sptr->get_size()), counts(
                comm_sptr->get_size())
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
        }
        grid_shape[0] = grid_size0;
        grid_shape[1] = grid_size1;
        grid_shape[2] = grid_size2;
        rectangular_grid_domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(physical_size, physical_offset,
                        grid_shape, true));
        decompose_1d(*comm_sptr, grid_shape[0], offsets, counts);
        int lower = offsets[comm_sptr->get_rank()];
        int upper = offsets[comm_sptr->get_rank()] + counts[comm_sptr->get_rank()];
        distributed_rectangular_grid_sptr = Distributed_rectangular_grid_sptr(
                new Distributed_rectangular_grid(rectangular_grid_domain_sptr,
                        lower, upper, comm_sptr));
    }

    ~Distributed_rectangular_grid_domain_fixture_periodic()
    {
    }

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Commxx_sptr comm_sptr;
    std::vector<int > offsets, counts;
    Rectangular_grid_domain_sptr rectangular_grid_domain_sptr;
    Distributed_rectangular_grid_sptr distributed_rectangular_grid_sptr;
};

struct Distributed_rectangular_grid_domain_fixture_rectangular
{
    Distributed_rectangular_grid_domain_fixture_rectangular() :
        physical_size(3), physical_offset(3), grid_shape(3), comm_sptr(
                new Commxx), offsets(comm_sptr->get_size()), counts(
                comm_sptr->get_size())
    {
        for (int i = 0; i < 3; ++i) {
            physical_offset[i] = domain_offset;
            physical_size[i] = domain_max - domain_min;
        }
        grid_shape[0] = grid_size0;
        grid_shape[1] = grid_size1;
        grid_shape[2] = grid_size2;
        rectangular_grid_domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(physical_size, physical_offset,
                        grid_shape, true));
        decompose_1d(*comm_sptr, grid_shape[0], offsets, counts);
        int lower = offsets[comm_sptr->get_rank()];
        int upper = offsets[comm_sptr->get_rank()] + counts[comm_sptr->get_rank()];
        std::string solver("rectangular");
        distributed_rectangular_grid_sptr = Distributed_rectangular_grid_sptr(
                new Distributed_rectangular_grid(rectangular_grid_domain_sptr,
                        lower, upper, comm_sptr, solver));
    }

    ~Distributed_rectangular_grid_domain_fixture_rectangular()
    {
    }

    std::vector<double > physical_size, physical_offset;
    std::vector<int > grid_shape;
    Commxx_sptr comm_sptr;
    std::vector<int > offsets, counts;
    Rectangular_grid_domain_sptr rectangular_grid_domain_sptr;
    Distributed_rectangular_grid_sptr distributed_rectangular_grid_sptr;
};



double
f(int i, int j, int k, int period)
{
    double ii = i;
    if (i < 0) {
        ii = i + period;
    }
    if (i >= period) {
        ii = i - period;
    }
    return 100.0 * (ii + 1) + (j + 1) + 0.01 * (k + 1);
}

void
zero(MArray3d_ref & a, std::vector<int > const & grid_shape, int lower_guard,
        int upper_guard)
{
    for (int i = lower_guard; i < upper_guard; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                a[i][j][k] = 0.0;
            }
        }
    }
}

void
fill(MArray3d_ref & a, std::vector<int > const & grid_shape, int lower,
        int upper, int period)
{
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                a[i][j][k] = f(i, j, k, period);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(construct, Distributed_rectangular_grid_domain_fixture)
{

}

BOOST_FIXTURE_TEST_CASE(construct_periodic, Distributed_rectangular_grid_domain_fixture_periodic)
{

}

BOOST_FIXTURE_TEST_CASE(get_upper_lower, Distributed_rectangular_grid_domain_fixture)
{
    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_lower(),
            offsets[comm_sptr->get_rank()]);
    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_upper(),
            offsets[comm_sptr->get_rank()] + counts[comm_sptr->get_rank()]);
}

BOOST_FIXTURE_TEST_CASE(fill_guards, Distributed_rectangular_grid_domain_fixture)
{
    int lower = distributed_rectangular_grid_sptr->get_lower();
    int upper = distributed_rectangular_grid_sptr->get_upper();
    int lower_guard = distributed_rectangular_grid_sptr->get_lower_guard();
    int upper_guard = distributed_rectangular_grid_sptr->get_upper_guard();

    zero(distributed_rectangular_grid_sptr->get_grid_points(), grid_shape,
            lower_guard, upper_guard);
    fill(distributed_rectangular_grid_sptr->get_grid_points(), grid_shape,
            lower, upper, grid_shape[0]);

    if (lower != lower_guard) {
        for (int i = lower_guard; i < lower; ++i) {
            for (int j = 0; j < grid_shape[1]; ++j) {
                for (int k = 0; k < grid_shape[2]; ++k) {
                    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                            0.0);
                }
            }
        }
    }
    if (upper != upper_guard) {
        for (int i = upper; i < upper_guard; ++i) {
            for (int j = 0; j < grid_shape[1]; ++j) {
                for (int k = 0; k < grid_shape[2]; ++k) {
                    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                            0.0);
                }
            }
        }
    }

    distributed_rectangular_grid_sptr->fill_guards();

    for (int i = lower_guard; i < upper_guard; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                BOOST_CHECK_CLOSE(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                        f(i,j,k, grid_shape[0]), tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(fill_guards_periodic,
        Distributed_rectangular_grid_domain_fixture_periodic)
{
    int lower = distributed_rectangular_grid_sptr->get_lower();
    int upper = distributed_rectangular_grid_sptr->get_upper();
    int lower_guard = distributed_rectangular_grid_sptr->get_lower_guard();
    int upper_guard = distributed_rectangular_grid_sptr->get_upper_guard();

    zero(distributed_rectangular_grid_sptr->get_grid_points(), grid_shape,
            lower_guard, upper_guard);
    fill(distributed_rectangular_grid_sptr->get_grid_points(), grid_shape,
            lower, upper, grid_shape[0]);
    if (lower != lower_guard) {
        for (int i = lower_guard; i < lower; ++i) {
            for (int j = 0; j < grid_shape[1]; ++j) {
                for (int k = 0; k < grid_shape[2]; ++k) {
                    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                            0.0);
                }
            }
        }
    }
    if (upper != upper_guard) {
        for (int i = upper; i < upper_guard; ++i) {
            for (int j = 0; j < grid_shape[1]; ++j) {
                for (int k = 0; k < grid_shape[2]; ++k) {
                    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                            0.0);
                }
            }
        }
    }

    distributed_rectangular_grid_sptr->fill_guards();
    for (int i = lower_guard; i < upper_guard; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                BOOST_CHECK_CLOSE(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                        f(i,j,k, grid_shape[0]), tolerance);
            }
        }
    }
}

BOOST_FIXTURE_TEST_CASE(fill_guards_empty_processors,
        Distributed_rectangular_grid_domain_with_empty_fixture)
{
    int lower = distributed_rectangular_grid_sptr->get_lower();
    int upper = distributed_rectangular_grid_sptr->get_upper();
    int lower_guard = distributed_rectangular_grid_sptr->get_lower_guard();
    int upper_guard = distributed_rectangular_grid_sptr->get_upper_guard();

    zero(distributed_rectangular_grid_sptr->get_grid_points(), grid_shape,
            lower_guard, upper_guard);
    fill(distributed_rectangular_grid_sptr->get_grid_points(), grid_shape,
            lower, upper, grid_shape[0]);

    if (lower != lower_guard) {
        for (int i = lower_guard; i < lower; ++i) {
            for (int j = 0; j < grid_shape[1]; ++j) {
                for (int k = 0; k < grid_shape[2]; ++k) {
                    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                            0.0);
                }
            }
        }
    }
    if (upper != upper_guard) {
        for (int i = upper; i < upper_guard; ++i) {
            for (int j = 0; j < grid_shape[1]; ++j) {
                for (int k = 0; k < grid_shape[2]; ++k) {
                    BOOST_CHECK_EQUAL(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                            0.0);
                }
            }
        }
    }

    distributed_rectangular_grid_sptr->fill_guards();

    for (int i = lower_guard; i < upper_guard; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                BOOST_CHECK_CLOSE(distributed_rectangular_grid_sptr->get_grid_points()[i][j][k],
                        f(i,j,k, grid_shape[0]), tolerance);
            }
        }
    }
}



BOOST_FIXTURE_TEST_CASE(construct_rectangular, Distributed_rectangular_grid_domain_fixture_rectangular)
{

}
