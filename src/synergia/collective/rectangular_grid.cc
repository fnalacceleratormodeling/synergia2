#include "rectangular_grid.h"

Rectangular_grid::Rectangular_grid(std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic_z, storage3d storage) :
    normalization(1.0), storage(storage)
{
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            physical_size, physical_offset, grid_shape, periodic_z));
    grid_points_sptr = boost::shared_ptr<Raw_MArray3d >(new Raw_MArray3d(
            boost::extents[grid_shape[0]][grid_shape[1]][grid_shape[2]],storage));
    grid_points_2dc_sptr = boost::shared_ptr<Raw_MArray2dc >(new Raw_MArray2dc(
            boost::extents[grid_shape[0]][grid_shape[1]]));
    grid_points_1d_sptr = boost::shared_ptr<Raw_MArray1d >(new Raw_MArray1d(
            boost::extents[grid_shape[2]]));
}

Rectangular_grid::Rectangular_grid(
        Rectangular_grid_domain_sptr rectangular_grid_domain_sptr, storage3d storage) :
    normalization(1.0), storage(storage)
{
    domain_sptr = rectangular_grid_domain_sptr;
    std::vector<int >
            grid_shape(rectangular_grid_domain_sptr->get_grid_shape());
    grid_points_sptr = boost::shared_ptr<Raw_MArray3d >(new Raw_MArray3d(
            boost::extents[grid_shape[0]][grid_shape[1]][grid_shape[2]],storage));
    grid_points_2dc_sptr = boost::shared_ptr<Raw_MArray2dc >(new Raw_MArray2dc(
            boost::extents[grid_shape[0]][grid_shape[1]]));
    grid_points_1d_sptr = boost::shared_ptr<Raw_MArray1d >(new Raw_MArray1d(
            boost::extents[grid_shape[2]]));
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
    return grid_points_sptr->m;
}

MArray3d_ref &
Rectangular_grid::get_grid_points()
{
    return grid_points_sptr->m;
}

MArray2dc_ref const&
Rectangular_grid::get_grid_points_2dc() const
{
    return grid_points_2dc_sptr->m;
}

MArray2dc_ref &
Rectangular_grid::get_grid_points_2dc()
{
    return grid_points_2dc_sptr->m;
}

MArray1d_ref const&
Rectangular_grid::get_grid_points_1d() const
{
    return grid_points_1d_sptr->m;
}

MArray1d_ref &
Rectangular_grid::get_grid_points_1d()
{
    return grid_points_1d_sptr->m;
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

storage3d
Rectangular_grid::get_storage() const
{
    return storage;
}

