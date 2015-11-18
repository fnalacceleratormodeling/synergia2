#include "rectangular_grid_domain.h"

Rectangular_grid_domain::Rectangular_grid_domain(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape, bool periodic_z) :
    physical_size(3), physical_offset(3), grid_shape(3), left(3), cell_size(3),
    periodic_z(periodic_z)
{
    this->physical_size = physical_size;
    this->physical_offset = physical_offset;
    this->grid_shape = grid_shape;

    for (int i = 0; i < 3; ++i) {
        left[i] = physical_offset[i] - physical_size[i] / 2.0;
        cell_size[i] = physical_size[i] / (1.0 * grid_shape[i]);
    }
}

Rectangular_grid_domain::Rectangular_grid_domain(
        std::vector<double > const & physical_size,
        std::vector<double > const & physical_offset,
        std::vector<int > const & grid_shape) :
    physical_size(3), physical_offset(3), grid_shape(3), left(3), cell_size(3),
    periodic_z(false)
{
    this->physical_size = physical_size;
    this->physical_offset = physical_offset;
    this->grid_shape = grid_shape;

    for (int i = 0; i < 3; ++i) {
        left[i] = physical_offset[i] - physical_size[i] / 2.0;
        cell_size[i] = physical_size[i] / (1.0 * grid_shape[i]);
    }
}

Rectangular_grid_domain::Rectangular_grid_domain(
        std::vector<double > const & physical_size,
        std::vector<int > const & grid_shape, bool periodic_z) :
    physical_size(3), physical_offset(3), grid_shape(3), left(3), cell_size(3),
    periodic_z(periodic_z)
{
    this->physical_size = physical_size;
    
    std::vector<double > offset(3,0.);
    this->physical_offset = offset;
    
    this->grid_shape = grid_shape;
    this->periodic_z = periodic_z;
   
    for (int i = 0; i < 3; ++i) {
        left[i] = offset[i] - physical_size[i] / 2.0;
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
    return periodic_z;
}
std::vector<double > const&
Rectangular_grid_domain::get_left() const
{
    return left;
}

void
Rectangular_grid_domain::get_cell_coordinates(int ix, int iy, int iz,
        double & x, double & y, double & z) const
{
    x = left[0] + cell_size[0] * (0.5 + ix);
    y = left[1] + cell_size[1] * (0.5 + iy);
    z = left[2] + cell_size[2] * (0.5 + iz);
}
