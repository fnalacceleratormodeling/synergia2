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
    bool in_domain = domain.get_leftmost_indices_offsets(z, y, x, iz, iy, ix,
            offz, offy, offx);
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

inline std::complex<double >
interpolate_rectangular_2d(double x, double y, double z,
        Rectangular_grid_domain const& domain, MArray2dc_ref const& a,
        MArray1d_ref const& b)
{
    // bi-linear interpolation
    int ix, iy, iz;
    double offx, offy, offz;
    bool in_domain = domain.get_leftmost_indices_offsets(x, y, z, ix, iy, iz,
            offx, offy, offz);
    std::complex<double > val;
    if (in_domain) {
        double line_density = ((1.0 - offz) * b[iz] + offz * b[iz + 1]);
        val = line_density * ((1.0 - offx) * (1.0 - offy) * a[ix][iy]
                + offx * (1.0 - offy) * a[ix + 1][iy]
                + (1.0 - offx) * offy * a[ix][iy + 1]
                + offx * offy * a[ix + 1][iy + 1]);
    }
    return val;
}

inline std::complex<double >
interpolate_rectangular_2d(Raw_MArray1d const& bin,
        std::vector<int > const& grid_shape, bool periodic_z,
        MArray2dc_ref const& a, MArray1d_ref const& b)
{
    // bi-linear interpolation
    int ix, iy, iz;
    double offx, offy, offz;
    double aoffx, aoffy, aoffz;
    ix = fast_int_floor(bin.m[0]);
    iy = fast_int_floor(bin.m[2]);
    iz = fast_int_floor(bin.m[4]);
    offx = bin.m[1];
    offy = bin.m[3];
    offz = bin.m[5];
    std::complex<double > val;
    if (((ix>=0) && (ix<grid_shape[0] - 1) && (iy>=0) && (iy<grid_shape[1] - 1))
        && (periodic_z || ((iz>=0) && (iz<grid_shape[2] - 1)))) {
        double line_density = ((1.0 - offz) * b[iz] + offz * b[iz + 1]);
        val = line_density * ((1.0 - offx) * (1.0 - offy) * a[ix][iy]
                + offx * (1.0 - offy) * a[ix + 1][iy]
                + (1.0 - offx) * offy * a[ix][iy + 1]
                + offx * offy * a[ix + 1][iy + 1]);
    }
    return val;
}

#endif /* INTERPOLATE_RECTANGULAR_ZYX_H_ */
