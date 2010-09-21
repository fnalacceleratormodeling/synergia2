#ifndef SPACE_CHARGE_3D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_3D_OPEN_HOCKNEY_H_
#include "components/simulation/operator.h"
#include "components/bunch/bunch.h"
#include "components/collective/rectangular_grid_domain.h"
#include "components/collective/rectangular_grid.h"
#include "components/collective/distributed_rectangular_grid.h"
#include "utils/commxx.h"
#include "utils/distributed_fft3d.h"

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_3d_open_hockney : public Collective_operator
{
private:
    std::vector<int > grid_shape, doubled_grid_shape;
    Rectangular_grid_domain_sptr domain_sptr, doubled_domain_sptr;
    bool periodic_z;
    double z_period;
    Distributed_fft3d_sptr distributed_fft3d_sptr;
    Commxx comm;
    double n_sigma;
public:
    Space_charge_3d_open_hockney(std::vector<int > const & grid_shape,
            bool periodic_z, Commxx const& comm, double z_period = 0.0,
            double n_sigma = 4.0);
    /// Note: Use Space_charge_3d_open_hockney::get_internal_grid_shape for
    /// Distributed_fft3d.
    Space_charge_3d_open_hockney(bool periodic_z,
            Distributed_fft3d_sptr distributed_fft3d_sptr,
            double z_period = 0.0, double n_sigma = 4.0);
    void
    update_domain(Bunch const& bunch);
    Rectangular_grid_domain_sptr
    get_domain_sptr();
    Rectangular_grid_domain_sptr
    get_doubled_domain_sptr();
    /// Returns local charge density on original grid in [C/m^3]
    Rectangular_grid_sptr
    get_local_charge_density(Bunch const& bunch);
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2(Rectangular_grid const& local_charge_density);
    /// Returns Green function on the doubled grid in [1/m^3]
    Distributed_rectangular_grid_sptr
    get_green_fn2();
    Distributed_rectangular_grid_sptr
    get_scalar_field2(Distributed_rectangular_grid & charge_density22,
            Distributed_rectangular_grid & green_fn2);
    Distributed_rectangular_grid_sptr
    extract_scalar_field(
            Distributed_rectangular_grid const& scalar_field2);
    virtual
    void
    apply(Bunch & bunch, Operators & step_operators);
    ~Space_charge_3d_open_hockney();
};

typedef boost::shared_ptr<Space_charge_3d_open_hockney >
        Space_charge_3d_open_hockney_sptr;

#endif /* SPACE_CHARGE_3D_OPEN_HOCKNEY_H_ */
