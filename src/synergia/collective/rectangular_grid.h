#ifndef RECTANGULAR_GRID_H_
#define RECTANGULAR_GRID_H_
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/utils/multi_array_typedefs.h"

class Rectangular_grid
{
private:
    Rectangular_grid_domain_sptr domain_sptr;
    boost::shared_ptr<MArray3d > grid_points_sptr;
    double normalization;
public:
    Rectangular_grid(std::vector<double > const & physical_size, std::vector<
            double > const & physical_offset,
            std::vector<int > const & grid_shape, bool periodic_z);
    Rectangular_grid(
            Rectangular_grid_domain_sptr rectangular_grid_domain_sptr);
    Rectangular_grid_domain_sptr
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr
    get_domain_sptr();
    MArray3d_ref const&
    get_grid_points() const;
    MArray3d_ref &
    get_grid_points();
    void
    set_normalization(double val);
    double
    get_normalization() const;
};

typedef boost::shared_ptr<Rectangular_grid> Rectangular_grid_sptr;

#endif /* RECTANGULAR_GRID_H_ */
