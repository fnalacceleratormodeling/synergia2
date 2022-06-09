#include "synergia/collective/deposit.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/distributed_fft3d.h"

#include <Kokkos_ScatterView.hpp>

namespace deposit_impl {

  KOKKOS_INLINE_FUNCTION
  int
  fast_int_floor_kokkos(const double x)
  {
    int ix = static_cast<int>(x);
    return x > 0.0 ? ix : ((x - ix == 0) ? ix : ix - 1);
  }

  KOKKOS_INLINE_FUNCTION
  void
  get_leftmost_indices_offset(double pos,
                              double left,
                              double inv_cell_size,
                              int& idx,
                              double& off)
  {
    double scaled_location = (pos - left) * inv_cell_size - 0.5;
    idx = fast_int_floor_kokkos(scaled_location);
    off = scaled_location - idx;
  }

  KOKKOS_INLINE_FUNCTION
  bool
  ingrid(int i, int g)
  {
    return i >= 0 && i < g;
  }

  KOKKOS_INLINE_FUNCTION
  bool
  ingrid(int ix, int iy, int iz, int gx, int gy, int gz)
  {
    // return ix>=0 && ix<gx && iy>=0 && iy<gy && iz>=0 && iz<gz;

    // exclude edges
    return ix > 0 && ix < gx - 1 && iy > 0 && iy < gy - 1 && iz > 0 &&
           iz < gz - 1;
  }

  struct rho_reducer {
    typedef double value_type[];

    const int value_count;
    ConstParticles p;
    karray2d_dev bin;
    int gx, gy, gz;
    double ihx, ihy, ihz;
    double lx, ly, lz;
    double w0;

    rho_reducer(ConstParticles const& p,
                karray2d_dev const& bin,
                std::array<int, 3> const& g,
                std::array<double, 3> const& h,
                std::array<double, 3> const& l,
                double w0)
      : value_count(g[0] * g[1] + g[2])
      , p(p)
      , bin(bin)
      , gx(g[0])
      , gy(g[1])
      , gz(g[2])
      , ihx(1.0 / h[0])
      , ihy(1.0 / h[1])
      , ihz(1.0 / h[2])
      , lx(l[0])
      , ly(l[1])
      , lz(l[2])
      , w0(w0)
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i, value_type sum) const
    {
      int ix, iy, iz;
      double offx, offy, offz;

      get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, offx);
      get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, offy);
      get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, offz);

      bin(i, 0) = ix;
      bin(i, 1) = offx;
      bin(i, 2) = iy;
      bin(i, 3) = offy;
      bin(i, 4) = iz;
      bin(i, 5) = offz;

      int cellz1 = iz;
      int cellz2 = cellz1 + 1;

      if (ingrid(cellz1, gz)) sum[gx * gy + cellz1] += (1.0 - offz) * ihz;
      if (ingrid(cellz2, gz)) sum[gx * gy + cellz2] += offz * ihz;

      if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) return;

      int cellx1, cellx2, celly1, celly2;
      cellx1 = ix;
      cellx2 = ix + 1;
      celly1 = iy;
      celly2 = iy + 1;

      double aoffx, aoffy;
      aoffx = 1. - offx;
      aoffy = 1. - offy;

      sum[cellx1 * gy + celly1] += w0 * aoffx * aoffy;
      sum[cellx1 * gx + celly2] += w0 * aoffx * offy;
      sum[cellx2 * gy + celly1] += w0 * offx * aoffy;
      sum[cellx2 * gx + celly2] += w0 * offx * offy;
    }
  };

  // move the rho data from double array to complex array
  struct rho_mover {
    karray1d_dev r0;
    karray1d_dev r1;
    int gx, gy;

    rho_mover(karray1d_dev const& rho_double,
              karray1d_dev const& rho_complex,
              std::array<int, 3> const& g)
      : r0(rho_double), r1(rho_complex), gx(g[0]), gy(g[1])
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      if (i < gx * gy) {
        r1(i * 2) = r0(i);
      } else {
        int z = i - gx * gy;
        r1(gx * gy * 2 + z) = r0(gx * gy + z);
      }
    }
  };

  // use atomic add
  struct atomic_rho_reducer {
    ConstParticles p;
    karray1d_atomic_dev rho;
    karray2d_dev bin;
    int gx, gy, gz;
    double ihx, ihy, ihz;
    double lx, ly, lz;
    double w0;

    atomic_rho_reducer(ConstParticles const& p,
                       karray1d_atomic_dev const& rho,
                       karray2d_dev const& bin,
                       std::array<int, 3> const& g,
                       std::array<double, 3> const& h,
                       std::array<double, 3> const& l,
                       double w0)
      : p(p)
      , rho(rho)
      , bin(bin)
      , gx(g[0])
      , gy(g[1])
      , gz(g[2])
      , ihx(1.0 / h[0])
      , ihy(1.0 / h[1])
      , ihz(1.0 / h[2])
      , lx(l[0])
      , ly(l[1])
      , lz(l[2])
      , w0(w0)
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      int ix, iy, iz;
      double offx, offy, offz;

      get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, offx);
      get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, offy);
      get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, offz);

      bin(i, 0) = ix;
      bin(i, 1) = offx;
      bin(i, 2) = iy;
      bin(i, 3) = offy;
      bin(i, 4) = iz;
      bin(i, 5) = offz;

      int cellz1 = iz;
      int cellz2 = cellz1 + 1;

      if (cellz1 >= 0 && cellz1 < gz)
        rho(gx * gy * 2 + cellz1) += (1.0 - offz) * ihz;
      if (cellz2 >= 0 && cellz2 < gz) rho(gx * gy * 2 + cellz2) += offz * ihz;

      if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) return;

      int cellx1, cellx2, celly1, celly2;
      cellx1 = ix;
      cellx2 = ix + 1;
      celly1 = iy;
      celly2 = iy + 1;

      double aoffx, aoffy;
      aoffx = 1. - offx;
      aoffy = 1. - offy;

      rho((cellx1 * gy + celly1) * 2) += w0 * aoffx * aoffy;
      rho((cellx1 * gx + celly2) * 2) += w0 * aoffx * offy;
      rho((cellx2 * gy + celly1) * 2) += w0 * offx * aoffy;
      rho((cellx2 * gx + celly2) * 2) += w0 * offx * offy;
    }
  };

  struct rho_zeroer {
    karray1d_dev rho;

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      rho(i) = 0.0;
    }
  };

  // use scatter view
  struct sv_rho_reducer {
    ConstParticles p;
    ConstParticleMasks masks;

    Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft> scatter;
    karray2d_dev bin;
    int gx, gy, gz;
    double ihx, ihy, ihz;
    double lx, ly, lz;
    double w0;

    sv_rho_reducer(
      ConstParticles const& p,
      ConstParticleMasks const& masks,
      Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft> const&
        scatter,
      karray2d_dev const& bin,
      std::array<int, 3> const& g,
      std::array<double, 3> const& h,
      std::array<double, 3> const& l,
      double w0)
      : p(p)
      , masks(masks)
      , scatter(scatter)
      , bin(bin)
      , gx(g[0])
      , gy(g[1])
      , gz(g[2])
      , ihx(1.0 / h[0])
      , ihy(1.0 / h[1])
      , ihz(1.0 / h[2])
      , lx(l[0])
      , ly(l[1])
      , lz(l[2])
      , w0(w0)
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      if (masks(i)) {
        auto access = scatter.access();

        int ix, iy, iz;
        double offx, offy, offz;

        get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, offx);
        get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, offy);
        get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, offz);

        bin(i, 0) = ix;
        bin(i, 1) = offx;
        bin(i, 2) = iy;
        bin(i, 3) = offy;
        bin(i, 4) = iz;
        bin(i, 5) = offz;

        int cellz1 = iz;
        int cellz2 = cellz1 + 1;

        if (cellz1 >= 0 && cellz1 < gz)
          access(gx * gy * 2 + cellz1) += (1.0 - offz) * ihz;
        if (cellz2 >= 0 && cellz2 < gz)
          access(gx * gy * 2 + cellz2) += offz * ihz;

        if (ix < 0 || ix > gx - 1 || iy < 0 || iy > gy - 1) return;

        int cellx1, cellx2, celly1, celly2;
        cellx1 = ix;
        cellx2 = ix + 1;
        celly1 = iy;
        celly2 = iy + 1;

        double aoffx, aoffy;
        aoffx = 1. - offx;
        aoffy = 1. - offy;

        access((cellx1 * gy + celly1) * 2) += w0 * aoffx * aoffy;
        access((cellx1 * gx + celly2) * 2) += w0 * aoffx * offy;
        access((cellx2 * gy + celly1) * 2) += w0 * offx * aoffy;
        access((cellx2 * gx + celly2) * 2) += w0 * offx * offy;
      }
    }
  };

  // use scatter view
  struct sv_zyx_rho_reducer_non_periodic {
    ConstParticles p;
    ConstParticleMasks masks;

    scatter_t scatter;

    int gx, gy, gz; // original grid size
    int dx, dy, dz; // dimensions of the grid
    double ihx, ihy, ihz;
    double lx, ly, lz;
    double w0;

    sv_zyx_rho_reducer_non_periodic(ConstParticles const& p,
                                    ConstParticleMasks const& masks,
                                    scatter_t const& scatter,
                                    std::array<int, 3> const& g,
                                    std::array<int, 3> const& d,
                                    std::array<double, 3> const& h,
                                    std::array<double, 3> const& l,
                                    double w0)
      : p(p)
      , masks(masks)
      , scatter(scatter)
      , gx(g[0])
      , gy(g[1])
      , gz(g[2])
      , dx(d[0])
      , dy(d[1])
      , dz(d[2])
      , ihx(1.0 / h[0])
      , ihy(1.0 / h[1])
      , ihz(1.0 / h[2])
      , lx(l[0])
      , ly(l[1])
      , lz(l[2])
      , w0(w0)
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      if (masks(i)) {
        auto access = scatter.access();

        int ix, iy, iz;
        double ox, oy, oz;

        get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, ox);
        get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, oy);
        get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, oz);

        double aox = 1.0 - ox;
        double aoy = 1.0 - oy;
        double aoz = 1.0 - oz;

        int base = iz * dx * dy;

        if (ingrid(ix, iy, iz, gx, gy, gz))
          access(base + iy * dx + ix) += w0 * aox * aoy * aoz;

        if (ingrid(ix + 1, iy, iz, gx, gy, gz))
          access(base + iy * dx + ix + 1) += w0 * ox * aoy * aoz;

        if (ingrid(ix, iy + 1, iz, gx, gy, gz))
          access(base + (iy + 1) * dx + ix) += w0 * aox * oy * aoz;

        if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz))
          access(base + (iy + 1) * dx + ix + 1) += w0 * ox * oy * aoz;

        base = (iz + 1) * dx * dy;

        if (ingrid(ix, iy, iz + 1, gx, gy, gz))
          access(base + iy * dx + ix) += w0 * aox * aoy * oz;

        if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz))
          access(base + iy * dx + ix + 1) += w0 * ox * aoy * oz;

        if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz))
          access(base + (iy + 1) * dx + ix) += w0 * aox * oy * oz;

        if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz))
          access(base + (iy + 1) * dx + ix + 1) += w0 * ox * oy * oz;
      }
    }
  };

  // use scatter view
  struct sv_xyz_rho_reducer_non_periodic {
    ConstParticles p;
    ConstParticleMasks masks;

    scatter_t scatter;

    int gx, gy, gz; // original grid size
    int dx, dy, dz; // dimensions of the grid
    double ihx, ihy, ihz;
    double lx, ly, lz;
    double w0;

    sv_xyz_rho_reducer_non_periodic(ConstParticles const& p,
                                    ConstParticleMasks const& masks,
                                    scatter_t const& scatter,
                                    std::array<int, 3> const& g,
                                    std::array<int, 3> const& d,
                                    std::array<double, 3> const& h,
                                    std::array<double, 3> const& l,
                                    double w0)
      : p(p)
      , masks(masks)
      , scatter(scatter)
      , gx(g[0])
      , gy(g[1])
      , gz(g[2])
      , dx(d[0])
      , dy(d[1])
      , dz(d[2])
      , ihx(1.0 / h[0])
      , ihy(1.0 / h[1])
      , ihz(1.0 / h[2])
      , lx(l[0])
      , ly(l[1])
      , lz(l[2])
      , w0(w0)
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      if (masks(i)) {
        auto access = scatter.access();

        int ix, iy, iz;
        double ox, oy, oz;

        get_leftmost_indices_offset(p(i, 0), lx, ihx, ix, ox);
        get_leftmost_indices_offset(p(i, 2), ly, ihy, iy, oy);
        get_leftmost_indices_offset(p(i, 4), lz, ihz, iz, oz);

        double aox = 1.0 - ox;
        double aoy = 1.0 - oy;
        double aoz = 1.0 - oz;

        int base = ix * dy * dz;

        if (ingrid(ix, iy, iz, gx, gy, gz))
          access(base + iy * dz + iz) += w0 * aox * aoy * aoz;

        if (ingrid(ix, iy, iz + 1, gx, gy, gz))
          access(base + iy * dz + iz + 1) += w0 * aox * aoy * oz;

        if (ingrid(ix, iy + 1, iz, gx, gy, gz))
          access(base + (iy + 1) * dz + iz) += w0 * aox * oy * aoz;

        if (ingrid(ix, iy + 1, iz + 1, gx, gy, gz))
          access(base + (iy + 1) * dz + iz + 1) += w0 * aox * oy * oz;

        base = (ix + 1) * dy * dz;

        if (ingrid(ix + 1, iy, iz, gx, gy, gz))
          access(base + iy * dz + iz) += w0 * ox * aoy * aoz;

        if (ingrid(ix + 1, iy, iz + 1, gx, gy, gz))
          access(base + iy * dz + iz + 1) += w0 * ox * aoy * oz;

        if (ingrid(ix + 1, iy + 1, iz, gx, gy, gz))
          access(base + (iy + 1) * dz + iz) += w0 * ox * oy * aoz;

        if (ingrid(ix + 1, iy + 1, iz + 1, gx, gy, gz))
          access(base + (iy + 1) * dz + iz + 1) += w0 * ox * oy * oz;
      }
    }
  };

}

karray1d_dev
deposit_charge_rectangular_2d_kokkos(Rectangular_grid_domain& domain,
                                     karray2d_dev& particle_bin,
                                     Bunch const& bunch)
{
  using deposit_impl::rho_mover;
  using deposit_impl::rho_reducer;

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  auto parts = bunch.get_local_particles();
  int nparts = bunch.get_local_num();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1]); // * h[2]);

  karray1d_dev rho_dbl("rho", g[0] * g[1] + g[2]);
  rho_reducer rr(parts, particle_bin, g, h, l, weight0);
  Kokkos::parallel_reduce(nparts, rr, rho_dbl);

  karray1d_dev rho_cplx("rho_cplx", g[0] * g[1] * 2 + g[2]);
  rho_mover rm(rho_dbl, rho_cplx, g);
  Kokkos::parallel_for(g[0] * g[1] + g[2], rm);
  Kokkos::fence();

  return rho_cplx;
}

karray1d_dev
deposit_charge_rectangular_2d_kokkos_atomic(Rectangular_grid_domain& domain,
                                            karray2d_dev& particle_bin,
                                            Bunch const& bunch)
{
  using deposit_impl::atomic_rho_reducer;

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  auto parts = bunch.get_local_particles();
  int nparts = bunch.get_local_num();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1]); // * h[2]);

  // double[x][y][2] + double[z]
  karray1d_dev rho_dev("rho", g[0] * g[1] * 2 + g[2]);
  karray1d_atomic_dev rho_atomic = rho_dev;

  atomic_rho_reducer rr(parts, rho_atomic, particle_bin, g, h, l, weight0);
  Kokkos::parallel_for(nparts, rr);
  Kokkos::fence();

  return rho_dev;
}

void
deposit_charge_rectangular_2d_kokkos_scatter_view(
  karray1d_dev& rho_dev,
  Rectangular_grid_domain& domain,
  karray2d_dev& particle_bin,
  Bunch const& bunch)
{
  using deposit_impl::rho_zeroer;
  using deposit_impl::sv_rho_reducer;

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  auto parts = bunch.get_local_particles();
  auto masks = bunch.get_local_particle_masks();
  int nparts = bunch.size();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1]); // * h[2]);

  // double[x][y][2] + double[z]
  if (rho_dev.extent(0) != g[0] * g[1] * 2 + g[2])
    throw std::runtime_error("wrong size for rho in deposit charge");

  // zero first
  rho_zeroer rz{rho_dev};
  Kokkos::parallel_for(rho_dev.extent(0), rz);
  Kokkos::fence();

  // deposit
  Kokkos::Experimental::ScatterView<double*, Kokkos::LayoutLeft> scatter(
    rho_dev);

  sv_rho_reducer rr(parts, masks, scatter, particle_bin, g, h, l, weight0);
  Kokkos::parallel_for(nparts, rr);
  Kokkos::Experimental::contribute(rho_dev, scatter);

  Kokkos::fence();
}

void
deposit_charge_rectangular_3d_kokkos_scatter_view(
  karray1d_dev& rho_dev,
  scatter_t& scatter_rho,
  Rectangular_grid_domain& domain,
  std::array<int, 3> const& dims,
  Bunch const& bunch)
{
  using namespace deposit_impl;

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  // g[0] (nx) needs to be padded to (2*(g[0]/2+1))
  // int padded_gx = Distributed_fft3d::get_padded_shape_real(g[0]);

  auto parts = bunch.get_local_particles();
  auto masks = bunch.get_local_particle_masks();
  int nparts = bunch.size();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1] * h[2]);

  if (rho_dev.extent(0) < g[0] * g[1] * g[2])
    throw std::runtime_error("insufficient size for rho in deposit charge");

  // zero first
  rho_zeroer rz{rho_dev};
  Kokkos::parallel_for(rho_dev.extent(0), rz);
  Kokkos::fence();

  scatter_rho.reset();
  // deposit
  sv_zyx_rho_reducer_non_periodic rr(
    parts, masks, scatter_rho, g, dims, h, l, weight0);

  Kokkos::parallel_for(nparts, rr);
  Kokkos::Experimental::contribute(rho_dev, scatter_rho);

  Kokkos::fence();
}

void
deposit_charge_rectangular_3d_kokkos_scatter_view_xyz(
  karray1d_dev& rho_dev,
  Rectangular_grid_domain& domain,
  std::array<int, 3> const& dims,
  Bunch const& bunch)
{
  using namespace deposit_impl;

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  // g[0] (nx) needs to be padded to (2*(g[0]/2+1))
  // int padded_gx = Distributed_fft3d::get_padded_shape_real(g[0]);

  auto parts = bunch.get_local_particles();
  auto masks = bunch.get_local_particle_masks();
  int nparts = bunch.size();

  double weight0 = (bunch.get_real_num() / bunch.get_total_num()) *
                   bunch.get_particle_charge() * pconstants::e /
                   (h[0] * h[1] * h[2]);

  if (rho_dev.extent(0) < g[0] * g[1] * g[2])
    throw std::runtime_error("insufficient size for rho in deposit charge");

  // zero first
  rho_zeroer rz{rho_dev};
  Kokkos::parallel_for(rho_dev.extent(0), rz);
  Kokkos::fence();

  // deposit
  scatter_t scatter(rho_dev);
  sv_xyz_rho_reducer_non_periodic rr(
    parts, masks, scatter, g, dims, h, l, weight0);

  Kokkos::parallel_for(nparts, rr);
  Kokkos::Experimental::contribute(rho_dev, scatter);

  Kokkos::fence();
}
