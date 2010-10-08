#include "rectangular_grid.h"

Rectangular_grid::Rectangular_grid(std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic_z) :
    normalization(1.0)
{
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape, periodic_z));
    grid_points_sptr = boost::shared_ptr<MArray3d >(new MArray3d(
            boost::extents[grid_shape[0]][grid_shape[1]][grid_shape[2]]));
}

Rectangular_grid::Rectangular_grid(
        Rectangular_grid_domain_sptr rectangular_grid_domain_sptr) :
    normalization(1.0)
{
    domain_sptr = rectangular_grid_domain_sptr;
    std::vector<int >
            grid_shape(rectangular_grid_domain_sptr->get_grid_shape());
    grid_points_sptr = boost::shared_ptr<MArray3d >(new MArray3d(
            boost::extents[grid_shape[0]][grid_shape[1]][grid_shape[2]]));
}

Rectangular_grid_domain_sptr
Rectangular_grid::get_domain_sptr() const
{
    return domain_sptr;
}

Rectangular_grid_domain_sptr
Rectangular_grid::get_domain_sptr()
{
    return domain_sptr;
}

MArray3d_ref const&
Rectangular_grid::get_grid_points() const
{
    return *grid_points_sptr;
}

MArray3d_ref &
Rectangular_grid::get_grid_points()
{
    return *grid_points_sptr;
}

void
Rectangular_grid::set_normalization(double val)
{
    normalization = val;
}

double
Rectangular_grid::get_normalization() const
{
    return normalization;
}
