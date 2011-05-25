#ifndef INTERPOLATE_RECTANGULAR_ZYX_H_
#define INTERPOLATE_RECTANGULAR_ZYX_H_

#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/distributed_rectangular_grid.h"

inline double
interpolate_rectangular_zyx(double x, double y, double z,
        Rectangular_grid_domain const& domain, MArray3d_ref const& a)
{
    // tri-linear interpolation
    int ix, iy, iz;
    double offx, offy, offz;
    bool in_domain = domain.get_leftmost_indices_offsets(z, y, x, iz, iy,
            ix, offz, offy, offx);
    double val;
    if (in_domain) {
        val = ((1.0 - offz) * (1.0 - offy) * (1.0 - offx) * a[iz][iy][ix]
                + offz * (1.0 - offy) * (1.0 - offx) * a[iz + 1][iy][ix] + (1.0
                - offz) * offy * (1.0 - offx) * a[iz][iy + 1][ix]
                + (1.0 - offz) * (1.0 - offy) * offx * a[iz][iy][ix + 1] + offz
                * offy * (1.0 - offx) * a[iz + 1][iy + 1][ix] + offz * (1.0
                - offy) * offx * a[iz + 1][iy][ix + 1] + (1.0 - offz) * offy
                * offx * a[iz][iy + 1][ix + 1] + offz * offy * offx
                * a[iz + 1][iy + 1][ix + 1]);
    } else {
        val = 0.0;
    }
    return val;
}

#endif /* INTERPOLATE_RECTANGULAR_ZYX_H_ */
