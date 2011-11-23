#ifndef DEPOSIT_H_
#define DEPOSIT_H_
#include "synergia/collective/rectangular_grid.h"
#include "synergia/bunch/bunch.h"

void
deposit_charge_rectangular_zyx(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);

void
deposit_charge_rectangular_xyz(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);


void
deposit_charge_rectangular_2d(Rectangular_grid & rho_grid, Bunch const& bunch,
        bool zero_first = true);


inline double
interpolate_rectangular_xyz(double x, double y, double z,
        Rectangular_grid_domain const& domain, MArray3d_ref const& a)
{
    // tri-linear interpolation
    int ix, iy, iz;
    double offx, offy, offz;
    bool in_domain = domain.get_leftmost_indices_offsets(x, y, z, ix, iy, iz,
            offx, offy, offz);
    double val;
    if (in_domain) {
        val = ((1.0 - offx) * (1.0 - offy) * (1.0 - offz) * a[ix][iy][iz]+
                  offx * (1.0 - offy) * (1.0 - offz) * a[ix + 1][iy][iz] + 
                  (1.0 - offx) * offy * (1.0 - offz) * a[ix][iy + 1][iz] +        
                + (1.0 - offx) * (1.0 - offy) * offz * a[ix][iy][iz + 1] + 
                       offx* offy * (1.0 - offz) * a[ix + 1][iy + 1][iz] + 
                       offx * (1.0- offy) * offz * a[ix + 1][iy][iz + 1] +
                       (1.0 - offx) * offy* offz * a[ix][iy + 1][iz + 1] +         
                           offx * offy * offz* a[ix + 1][iy + 1][iz + 1]);
    
    } else {
        val = 0.0;
    }
    return val;    
}        
        
#endif /* DEPOSIT_H_ */
