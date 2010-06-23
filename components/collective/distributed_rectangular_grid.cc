#include "distributed_rectangular_grid.h"

Distributed_rectangular_grid::Distributed_rectangular_grid(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic_z, int lower_z,
        int upper_z)
{
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape, periodic_z));
    grid_points_sptr
            = boost::shared_ptr<MArray3d >(
                    new MArray3d(
                            boost::extents[extent_range(lower_z, upper_z)][grid_shape[1]][grid_shape[2]]));
}

Distributed_rectangular_grid::Distributed_rectangular_grid(
        Rectangular_grid_domain_sptr const& Rectangular_grid_domain_sptr,
        int lower_z, int upper_z)
{
    domain_sptr = Rectangular_grid_domain_sptr;
    std::vector<int >
            grid_shape(Rectangular_grid_domain_sptr->get_grid_shape());
    grid_points_sptr
            = boost::shared_ptr<MArray3d >(
                    new MArray3d(
                            boost::extents[extent_range(lower_z, upper_z)][grid_shape[1]][grid_shape[2]]));
}

Rectangular_grid_domain_sptr &
Distributed_rectangular_grid::get_domain_sptr()
{
    return domain_sptr;
}

MArray3d_ref &
Distributed_rectangular_grid::get_grid_points()
{
    return *grid_points_sptr;
}
