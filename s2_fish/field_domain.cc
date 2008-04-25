#include "field_domain.h"
#include "array_nd/vector_helper.h"
#include <cmath>

void 
Field_domain::construct(const std::vector<double> &physical_size,
    const std::vector<double> &physical_offset,
    const std::vector<int> &grid_shape,
    const std::vector<bool> &periodic)
{
    this->physical_size = physical_size;
    this->physical_offset = physical_offset;
    this->grid_shape = grid_shape;
    this->periodic = periodic;
    
    left.resize(3);
    cell_size.resize(3);
    for (int i = 0; i < 3; ++i) {
        left[i] = physical_offset[i] - physical_size[i] / 2.0;
        cell_size[i] = physical_size[i] / (grid_shape[i] - 1.0);
    }
}

Field_domain::Field_domain()
{
    std::vector<double> zero_d = vector3<double>(0.0,0.0,0.0);
    std::vector<int> two_i = vector3<int>(2,2,2); 
    std::vector<bool> zero_b = vector3<bool>(false,false,false);
    construct(zero_d,zero_d,two_i,zero_b);
}

Field_domain::Field_domain(const std::vector<double> &physical_size,
    const std::vector<double> &physical_offset,
    const std::vector<int> &grid_shape,
    const std::vector<bool> &periodic)
{
    construct(physical_size,physical_offset,grid_shape,periodic);
}

void
Field_domain::set_params(const std::vector<double> &physical_size,
    const std::vector<double> &physical_offset,
    const std::vector<int> &grid_shape,
    const std::vector<bool> &periodic)
{
    construct(physical_size,physical_offset,grid_shape,periodic);
}

std::vector<int>
Field_domain::get_grid_shape() const
{
    return grid_shape;
}

std::vector<double>
Field_domain::get_physical_size() const
{
    return physical_size;
}

std::vector<double>
Field_domain::get_cell_size() const
{
    return cell_size;
}

std::vector<bool>
Field_domain::get_periodic() const
{
    return periodic;
}

void
Field_domain::get_leftmost_indices_offsets(double c0, double c1, double c2,
        std::vector<int> &indices, std::vector<double> &offsets) const
{
    double scaled_location;
    
    scaled_location = (c0 - left[0])/cell_size[0];
    offsets[0] = scaled_location - static_cast<int>(floor(scaled_location));
    indices[0] = static_cast<int>(floor(scaled_location));

    scaled_location = (c1 - left[1])/cell_size[1];
    offsets[1] = scaled_location - static_cast<int>(floor(scaled_location));
    indices[1] = static_cast<int>(floor(scaled_location));

    scaled_location = (c2 - left[2])/cell_size[2];
    offsets[2] = scaled_location - static_cast<int>(floor(scaled_location));
    indices[2] = static_cast<int>(floor(scaled_location));
}
