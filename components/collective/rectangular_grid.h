#ifndef RECTANGULAR_GRID_H_
#define RECTANGULAR_GRID_H_
#include "components/collective/rectangular_grid_domain.h"
#include "utils/multi_array_typedefs.h"

class Rectangular_grid
{
private:
    Rectangular_grid_domain_sptr domain_sptr;
    boost::shared_ptr<MArray3d > grid_points_sptr;
public:
    Rectangular_grid(std::vector<double > const & physical_size, std::vector<
            double > const & physical_offset,
            std::vector<int > const & grid_shape, bool periodic_z);
    Rectangular_grid(
            Rectangular_grid_domain_sptr const& rectangular_grid_domain_sptr);
    Rectangular_grid_domain_sptr &
    get_domain_sptr();
    MArray3d_ref &
    get_grid_points();
};

typedef boost::shared_ptr<Rectangular_grid> Rectangular_grid_sptr;

#endif /* RECTANGULAR_GRID_H_ */
