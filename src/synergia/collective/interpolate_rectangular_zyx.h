#ifndef INTERPOLATE_RECTANGULAR_ZYX_H_
#define INTERPOLATE_RECTANGULAR_ZYX_H_

#include "synergia/collective/rectangular_grid.h"

inline double
interpolate_rectangular_zyx(double x, double y, double z,
        Rectangular_grid const& f)
{
    // tri-linear interpolation
    int ix, iy, iz;
    double offx, offy, offz;
    f.get_domain_sptr()->get_leftmost_indices_offsets(z, y, x, iz, iy, ix,
            offz, offy, offx);
    MArray3d_ref a(f.get_grid_points());
    double val;
    if ((iz < 0) || (iz >= f.get_domain_sptr()->get_grid_shape()[0] - 1) || (iy
            < 0) || (iy >= f.get_domain_sptr()->get_grid_shape()[1] - 1) || (ix
            < 0) || (ix >= f.get_domain_sptr()->get_grid_shape()[2] - 1)) {
        val = 0.0;
    } else {
        val = ((1.0 - offz) * (1.0 - offy) * (1.0 - offx) * a[iz][iy][ix]
                + offz * (1.0 - offy) * (1.0 - offx) * a[iz + 1][iy][ix] + (1.0
                - offz) * offy * (1.0 - offx) * a[iz][iy + 1][ix]
                + (1.0 - offz) * (1.0 - offy) * offx * a[iz][iy][ix + 1] + offz
                * offy * (1.0 - offx) * a[iz + 1][iy + 1][ix] + offz * (1.0
                - offy) * offx * a[iz + 1][iy][ix + 1] + (1.0 - offz) * offy
                * offx * a[iz][iy + 1][ix + 1] + offz * offy * offx
                * a[iz + 1][iy + 1][ix + 1]);
    }
    return val;
}

#endif /* INTERPOLATE_RECTANGULAR_ZYX_H_ */
