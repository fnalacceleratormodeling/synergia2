#ifndef DEPOSIT_H_
#define DEPOSIT_H_
#include "synergia/collective/rectangular_grid.h"
#include "synergia/bunch/bunch.h"

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

        
#endif /* DEPOSIT_H_ */
