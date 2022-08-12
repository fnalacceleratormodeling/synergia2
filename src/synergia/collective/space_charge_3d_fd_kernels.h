#ifndef SPACE_CHARGE_3D_FD_KERNELS_H_
#define SPACE_CHARGE_3D_FD_KERNELS_H_

#include "utils.h"

#include "synergia/bunch/bunch.h"
#include "synergia/utils/kokkos_views.h"

struct alg_force_extractor {
  karray1d_dev phi2;
  karray1d_dev enx;
  karray1d_dev eny;
  karray1d_dev enz;

  int gx, gy, gz;
  double ihx, ihy, ihz;
  double igxgy;
  double igx;

  alg_force_extractor(karray1d_dev const& phi2,
                      karray1d_dev const& enx,
                      karray1d_dev const& eny,
                      karray1d_dev const& enz,
                      std::array<int, 3> const& g,
                      std::array<double, 3> const& h)
    : phi2(phi2)
    , enx(enx)
    , eny(eny)
    , enz(enz)
    , gx(g[0])
    , gy(g[1])
    , gz(g[2])
    , ihx(0.5 / h[0])
    , ihy(0.5 / h[1])
    , ihz(0.5 / h[2])
    , igxgy(1.0 / (gx * gy))
    , igx(1.0 / gx)
  {}

  KOKKOS_INLINE_FUNCTION
  void
  operator()(const int i) const
  {
    int iz = i * igxgy;
    int iy = (i - iz * gx * gy) * igx;
    int ix = i - iz * gx * gy - iy * gx;

    int ixl, ixr, iyl, iyr, izl, izr;

    double idx, idy, idz;

    // all boundaries will be skipped (set to 0)
    if (ix == 0 || ix == gx - 1 || iy == 0 || iy == gy - 1 || iz == 0 ||
        iz == gz - 1)
      return;

    ixl = ix - 1;
    ixr = ix + 1;
    iyl = iy - 1;
    iyr = iy + 1;
    izl = iz - 1;
    izr = iz + 1;

    idx = ihx;
    idy = ihy;
    idz = ihz;

    int idx_r = iz * gx * gy + iy * gx + ixr;
    int idx_l = iz * gx * gy + iy * gx + ixl;
    enx(i) = -(phi2(idx_r) - phi2(idx_l)) * idx;

    idx_r = iz * gx * gy + iyr * gx + ix;
    idx_l = iz * gx * gy + iyl * gx + ix;
    eny(i) = -(phi2(idx_r) - phi2(idx_l)) * idy;

    idx_r = izr * gx * gy + iy * gx + ix;
    idx_l = izl * gx * gy + iy * gx + ix;
    auto val_right = phi2(idx_r) * 1.0;
    auto val_left = phi2(idx_l) * 1.0;
    enz(i) = -(phi2(idx_r) - phi2(idx_l)) * idz;
  }
};

struct alg_kicker {
  Particles parts;
  ConstParticleMasks masks;

  karray1d_dev enx;
  karray1d_dev eny;
  karray1d_dev enz;

  int gx, gy, gz;
  double ihx, ihy, ihz;
  double lx, ly, lz;
  double factor, pref, m;

  alg_kicker(Particles parts,
             ConstParticleMasks masks,
             karray1d_dev const& enx,
             karray1d_dev const& eny,
             karray1d_dev const& enz,
             std::array<int, 3> const& g,
             std::array<double, 3> const& h,
             std::array<double, 3> const& l,
             double factor,
             double pref,
             double m)
    : parts(parts)
    , masks(masks)
    , enx(enx)
    , eny(eny)
    , enz(enz)
    , gx(g[0])
    , gy(g[1])
    , gz(g[2])
    , ihx(1.0 / h[0])
    , ihy(1.0 / h[1])
    , ihz(1.0 / h[2])
    , lx(l[0])
    , ly(l[1])
    , lz(l[2])
    , factor(factor)
    , pref(pref)
    , m(m)
  {}

  KOKKOS_INLINE_FUNCTION
  void
  operator()(const int i) const
  {
    if (masks(i)) {
      int ix, iy, iz;
      double ox, oy, oz;

      get_leftmost_indices_offset(parts(i, 0), lx, ihx, ix, ox);
      get_leftmost_indices_offset(parts(i, 2), ly, ihy, iy, oy);
      get_leftmost_indices_offset(parts(i, 4), lz, ihz, iz, oz);

      double aox = 1.0 - ox;
      double aoy = 1.0 - oy;
      double aoz = 1.0 - oz;

      if ((ix >= 0 && ix < gx - 1) && (iy >= 0 && iy < gy - 1) &&
          ((iz >= 0 && iz < gz - 1) /* || periodic_z */)) {
        double val = 0;
        int base = 0;

        // enz
        base = iz * gx * gy + iy * gx + ix;
        val = aoz * aoy * aox * enz(base);         // z, y, x
        val += aoz * aoy * ox * enz(base + 1);     // z, y, x+1
        val += aoz * oy * aox * enz(base + gx);    // z, y+1, x
        val += aoz * oy * ox * enz(base + gx + 1); // z, y+1, x+1

        base = (iz + 1) * gx * gy + iy * gx + ix;
        val += oz * aoy * aox * enz(base);        // z+1, y, x
        val += oz * aoy * ox * enz(base + 1);     // z+1, y, x+1
        val += oz * oy * aox * enz(base + gx);    // z+1, y+1, x
        val += oz * oy * ox * enz(base + gx + 1); // z+1, y+1, x+1

        double p = pref + parts(i, 5) * pref;
        double Eoc_i = std::sqrt(p * p + m * m);
        double Eoc_f = Eoc_i + factor * (-pref) * val;
        double dpop = (std::sqrt(Eoc_f * Eoc_f - m * m) -
                       std::sqrt(Eoc_i * Eoc_i - m * m)) /
                      pref;

        parts(i, 5) += dpop;

        // eny
        base = iz * gx * gy + iy * gx + ix;
        val = aoz * aoy * aox * eny(base);         // z, y, x
        val += aoz * aoy * ox * eny(base + 1);     // z, y, x+1
        val += aoz * oy * aox * eny(base + gx);    // z, y+1, x
        val += aoz * oy * ox * eny(base + gx + 1); // z, y+1, x+1

        base = (iz + 1) * gx * gy + iy * gx + ix;
        val += oz * aoy * aox * eny(base);        // z+1, y, x
        val += oz * aoy * ox * eny(base + 1);     // z+1, y, x+1
        val += oz * oy * aox * eny(base + gx);    // z+1, y+1, x
        val += oz * oy * ox * eny(base + gx + 1); // z+1, y+1, x+1

        parts(i, 3) += factor * val;

        // enx
        base = iz * gx * gy + iy * gx + ix;
        val = aoz * aoy * aox * enx(base);         // z, y, x
        val += aoz * aoy * ox * enx(base + 1);     // z, y, x+1
        val += aoz * oy * aox * enx(base + gx);    // z, y+1, x
        val += aoz * oy * ox * enx(base + gx + 1); // z, y+1, x+1

        base = (iz + 1) * gx * gy + iy * gx + ix;
        val += oz * aoy * aox * enx(base);        // z+1, y, x
        val += oz * aoy * ox * enx(base + 1);     // z+1, y, x+1
        val += oz * oy * aox * enx(base + gx);    // z+1, y+1, x
        val += oz * oy * ox * enx(base + gx + 1); // z+1, y+1, x+1

        parts(i, 1) += factor * val;
      }
    }
  }
};

#endif
