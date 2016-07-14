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

inline double
interpolate_rectangular_xyz(double x, double y, double z,
        Rectangular_grid_domain const& domain, MArray3d_ref const& a)
{
    // tri-linear interpolation
    int ix, iy, iz, iz1, izp1;
    double offx, offy, offz;
    bool in_domain = domain.get_leftmost_indices_offsets(x, y, z, ix, iy, iz1,
            offx, offy, offz);
    double val;
    if (in_domain) {
        int period= domain.get_grid_shape()[2];       
        iz=(iz1%period)*int(iz1>=0)+  (period - 1 - ((-iz1 - 1) % period))*int(iz1 < 0); 
        //izp1=((iz1+1)%period)*int((iz1+1)>=0)+  (period - 1 - ((-iz1-2) % period))*int((iz1+1) < 0); 
        izp1=(iz+1)%period;
        val = ((1.0 - offx) * (1.0 - offy) * (1.0 - offz) * a[ix][iy][iz]+
                  offx * (1.0 - offy) * (1.0 - offz) * a[ix + 1][iy][iz] + 
                  (1.0 - offx) * offy * (1.0 - offz) * a[ix][iy + 1][iz] +        
                + (1.0 - offx) * (1.0 - offy) * offz * a[ix][iy][izp1] + 
                       offx* offy * (1.0 - offz) * a[ix + 1][iy + 1][iz] + 
                       offx * (1.0- offy) * offz * a[ix + 1][iy][izp1] +
                       (1.0 - offx) * offy* offz * a[ix][iy + 1][izp1] +         
                           offx * offy * offz* a[ix + 1][iy + 1][izp1]);
    
    } else {
        val = 0.0;
    }
    return val;    
}        
        
#endif /* DEPOSIT_H_ */
