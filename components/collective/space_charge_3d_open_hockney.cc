#include "space_charge_3d_open_hockney.h"
#include "components/bunch/diagnostics.h"
#include "deposit.h"

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        std::vector<int > const & grid_shape, bool periodic_z,
        Commxx const& comm, double n_sigma) :
    Collective_operator("space charge"), grid_shape(3), periodic_z(periodic_z),
            comm(comm), n_sigma(n_sigma)
{
    this->grid_shape = get_internal_grid_shape(grid_shape);
    distributed_fft3d_sptr = Distributed_fft3d_sptr(new Distributed_fft3d(
            grid_shape, comm));
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(bool periodic_z,
        Distributed_fft3d_sptr const& distributed_fft3d_sptr, double n_sigma) :
    Collective_operator("space charge"), grid_shape(3), periodic_z(periodic_z),
            distributed_fft3d_sptr(distributed_fft3d_sptr), comm(
                    distributed_fft3d_sptr->get_comm()), n_sigma(n_sigma)
{
    grid_shape = distributed_fft3d_sptr->get_shape();
}

std::vector<int >
Space_charge_3d_open_hockney::get_internal_grid_shape(
        std::vector<int > const& external_grid_shape)
{
    std::vector<int > internal_grid_shape(3);

    internal_grid_shape[0] = 2 * external_grid_shape[2];
    internal_grid_shape[1] = 2 * external_grid_shape[1];
    internal_grid_shape[2] = 2 * external_grid_shape[0];

    return internal_grid_shape;
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_domain_sptr(Bunch const& bunch)
{
    Diagnostics diagnostics(bunch);
    std::vector<double > size(3);
    std::vector<double > offset(3);
    offset[0] = diagnostics.get_mean()[Bunch::z];
    size[0] = n_sigma * diagnostics.get_std()[Bunch::z];
    offset[1] = diagnostics.get_mean()[Bunch::y];
    size[1] = n_sigma * diagnostics.get_std()[Bunch::y];
    offset[0] = diagnostics.get_mean()[Bunch::z];
    size[0] = n_sigma * diagnostics.get_std()[Bunch::z];
    return Rectangular_grid_domain_sptr(new Rectangular_grid_domain(size,
            offset, grid_shape, periodic_z));
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    Rectangular_grid_sptr local_rho_sptr(new Rectangular_grid(get_domain_sptr(
            bunch)));
    deposit_charge_rectangular(*local_rho_sptr, bunch);
    return local_rho_sptr;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_charge_density(
        Rectangular_grid_sptr & local_charge_density_sptr)
{

}

void
Space_charge_3d_open_hockney::apply(Bunch & bunch, Operators & step_operators)
{
}

Space_charge_3d_open_hockney::~Space_charge_3d_open_hockney()
{

}

