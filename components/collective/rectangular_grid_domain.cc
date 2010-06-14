#include "rectangular_grid_domain.h"
#include <cmath>

Rectangular_grid_domain::Rectangular_grid_domain(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic_z) :
    physical_size(3), physical_offset(3), grid_shape(3), left(3), cell_size(3)
{
    this->physical_size = physical_size;
    this->physical_offset = physical_offset;
    this->grid_shape = grid_shape;
    this->periodic = periodic;

    for (int i = 0; i < 3; ++i) {
        left[i] = physical_offset[i] - physical_size[i] / 2.0;
        cell_size[i] = physical_size[i] / (1.0 * grid_shape[i]);
    }
}

std::vector<double > const&
Rectangular_grid_domain::get_physical_size() const
{
    return physical_size;
}

std::vector<double > const&
Rectangular_grid_domain::get_physical_offset() const
{
    return physical_offset;
}

std::vector<int > const&
Rectangular_grid_domain::get_grid_shape() const
{
    return grid_shape;
}

std::vector<double > const&
Rectangular_grid_domain::get_cell_size() const
{
    return cell_size;
}

bool
Rectangular_grid_domain::is_periodic() const
{
    return periodic;
}

void
Rectangular_grid_domain::get_leftmost_indices_offsets(double x, double y,
        double z, int & ix, int & iy, int & iz, double & offx, double & offy,
        double & offz) const
{
    double scaled_location;

    scaled_location = (x - left[0])/cell_size[0];
    ix = static_cast<int>(floor(scaled_location));
    offx = scaled_location - ix;

    scaled_location = (y - left[1])/cell_size[1];
    iy = static_cast<int>(floor(scaled_location));
    offy = scaled_location - iy;

    scaled_location = (z - left[2])/cell_size[2];
    iz = static_cast<int>(floor(scaled_location));
    offz = scaled_location - iz;
}
