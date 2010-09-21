#include "distributed_rectangular_grid.h"

void
Distributed_rectangular_grid::construct(int lower, int upper)
{
    std::vector<int > grid_shape(domain_sptr->get_grid_shape());
    this->lower = lower;
    this->upper = upper;
    if ((lower == 0) && (!domain_sptr->is_periodic())) {
        lower_guard = 0;
    } else {
        lower_guard = lower - 1;
    }
    if ((upper == grid_shape[0]) && (!domain_sptr->is_periodic())) {
        upper_guard = grid_shape[0];
    } else {
        upper_guard = upper + 1;
    }
    grid_points_sptr
            = boost::shared_ptr<MArray3d >(
                    new MArray3d(boost::extents[extent_range(lower_guard,
                            upper_guard)][grid_shape[1]][grid_shape[2]]));
    normalization = 1.0;
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic, int lower,
        int upper)
{
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape, periodic));
    construct(lower, upper);
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        Rectangular_grid_domain_sptr const& rectangular_grid_domain_sptr,
        int lower, int upper)
{
    domain_sptr = rectangular_grid_domain_sptr;
    construct(lower, upper);
}

Rectangular_grid_domain_sptr const&
Distributed_rectangular_grid::get_domain_sptr() const
{
    return domain_sptr;
}

Rectangular_grid_domain_sptr &
Distributed_rectangular_grid::get_domain_sptr()
{
    return domain_sptr;
}

int
Distributed_rectangular_grid::get_lower() const
{
    return lower;
}

int
Distributed_rectangular_grid::get_upper() const
{
    return upper;
}

int
Distributed_rectangular_grid::get_lower_guard() const
{
    return lower_guard;
}

int
Distributed_rectangular_grid::get_upper_guard() const
{
    return upper_guard;
}

MArray3d_ref const&
Distributed_rectangular_grid::get_grid_points() const
{
    return *grid_points_sptr;
}

MArray3d_ref &
Distributed_rectangular_grid::get_grid_points()
{
    return *grid_points_sptr;
}

void
Distributed_rectangular_grid::set_normalization(double val)
{
    normalization = val;
}

double
Distributed_rectangular_grid::get_normalization() const
{
    return normalization;
}
