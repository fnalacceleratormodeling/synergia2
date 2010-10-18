#ifndef RECTANGULAR_GRID_DOMAIN_H_
#define RECTANGULAR_GRID_DOMAIN_H_

#include <vector>
#include <boost/shared_ptr.hpp>

class Rectangular_grid_domain
{
private:
    std::vector<double > physical_size;
    std::vector<double > physical_offset;
    std::vector<int > grid_shape;
    bool periodic_z;

    std::vector<double > left;
    std::vector<double > cell_size;
public:
    Rectangular_grid_domain(std::vector<double > const & physical_size,
            std::vector<double > const & physical_offset,
            std::vector<int > const & grid_shape, bool periodic_z);

    std::vector<double > const&
    get_physical_size() const;
    std::vector<double > const&
    get_physical_offset() const;
    std::vector<int > const&
    get_grid_shape() const;
    std::vector<double > const&
    get_cell_size() const;
    bool
    is_periodic() const;
    void
    get_leftmost_indices_offsets(double x, double y, double z, int &ix,
            int &iy, int&iz, double &offx, double &offy, double &offz) const;
};

typedef boost::shared_ptr<Rectangular_grid_domain > Rectangular_grid_domain_sptr;

#endif /* RECTANGULAR_GRID_DOMAIN_H_ */
