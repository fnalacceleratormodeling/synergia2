#ifndef DISTRIBUTED_RECTANGULAR_GRID_H_
#define DISTRIBUTED_RECTANGULAR_GRID_H_
#include "components/collective/rectangular_grid_domain.h"
#include "utils/multi_array_typedefs.h"

class Distributed_rectangular_grid
{
private:
    Rectangular_grid_domain_sptr domain_sptr;
    boost::shared_ptr<MArray3d > grid_points_sptr;
    bool periodic;
    int lower, upper;
    int lower_guard, upper_guard;
    double normalization;
    void
    construct(int lower, int upper);
public:
    Distributed_rectangular_grid(std::vector<double > const & physical_size,
            std::vector<double > const & physical_offset,
            std::vector<int > const & grid_shape, bool periodic, int lower,
            int upper);
    Distributed_rectangular_grid(
            Rectangular_grid_domain_sptr const& rectangular_grid_domain_sptr,
            int lower, int upper);
    Rectangular_grid_domain_sptr const&
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr &
    get_domain_sptr();
    bool
    get_periodic() const;
    int
    get_lower() const;
    int
    get_upper() const;
    int
    get_lower_guard() const;
    int
    get_upper_guard() const;
    MArray3d_ref const&
    get_grid_points() const;
    MArray3d_ref &
    get_grid_points();
    void
    set_normalization(double val);
    double
    get_normalization() const;
};

typedef boost::shared_ptr<Distributed_rectangular_grid > Distributed_rectangular_grid_sptr;
#endif /* DISTRIBUTED_RECTANGULAR_GRID_H_ */
