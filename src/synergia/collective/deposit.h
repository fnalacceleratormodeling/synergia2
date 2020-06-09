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

void
deposit_charge_rectangular_2d_kokkos_scatter_view(
        karray1d_dev & rho_dev,
        Rectangular_grid_domain & domain,
        karray2d_dev & particle_bin, 
        Bunch const & bunch );

void
deposit_charge_rectangular_3d_kokkos_scatter_view(
        karray1d_dev& rho_dev,
        Rectangular_grid_domain& domain,
        std::array<int, 3> const& dims,
        Bunch const& bunch );

#ifdef Kokkos_ENABLE_OPENMP
void
deposit_charge_rectangular_2d_omp_reduce( 
        karray1d_dev & rho_dev,
        Rectangular_grid_domain & domain,
        karray2d_dev & bin, 
        Bunch const & bunch );

void
deposit_charge_rectangular_3d_omp_reduce(
        karray1d_dev& rho_dev,
        Rectangular_grid_domain& domain,
        std::array<int, 3> const& dims,
        Bunch const& bunch );
#endif




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
