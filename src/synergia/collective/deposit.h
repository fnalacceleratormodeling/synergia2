#ifndef DEPOSIT_H_
#define DEPOSIT_H_

#include "synergia/bunch/bunch.h"
#include "synergia/collective/rectangular_grid_domain.h"

karray1d_dev
deposit_charge_rectangular_2d_kokkos(
        Rectangular_grid_domain & domain,
        karray2d_dev & particle_bin, 
        Bunch const & bunch );

karray1d_dev
deposit_charge_rectangular_2d_kokkos_atomic(
        Rectangular_grid_domain & domain,
        karray2d_dev & particle_bin, 
        Bunch const & bunch );



#if 0
void
deposit_charge_rectangular_zyx(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);

void
deposit_charge_rectangular_zyx_omp_reduce(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);

void
deposit_charge_rectangular_zyx_omp_interleaved(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);

void
deposit_charge_rectangular_xyz(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);


void
deposit_charge_rectangular_2d(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);

void
deposit_charge_rectangular_2d(Rectangular_grid & rho_grid,
        Raw_MArray2d & particle_bin, Bunch const& bunch, bool zero_first = true);

void
deposit_charge_rectangular_2d_omp_reduce(Rectangular_grid & rho_grid,
        Raw_MArray2d & particle_bin, Bunch const& bunch, bool zero_first = true);
#endif

        
#endif /* DEPOSIT_H_ */
