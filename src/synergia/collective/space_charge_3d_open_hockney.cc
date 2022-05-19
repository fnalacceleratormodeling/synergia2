
#include "space_charge_3d_open_hockney.h"
#include "deposit.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"

#include "synergia/utils/simple_timer.h"

using mconstants::pi;

namespace {
  void
  print_grid(Logger& logger,
             karray1d_dev const& grid,
             int x0,
             int x1,
             int y0,
             int y1,
             int z0,
             int z1,
             int gx,
             int gy,
             int gz,
             int off = 0)
  {
    karray1d_hst hgrid = Kokkos::create_mirror_view(grid);
    Kokkos::deep_copy(hgrid, grid);

    double sum = 0;

    int dim = grid.extent(0);
    for (int i = 0; i < dim; ++i) { sum += fabs(hgrid(i)); }

#if 0
        for(int x=0; x<gx; ++x)
            for(int y=0; y<gy; ++y)
                sum += hgrid((x*gy + y)*2 + off);
#endif

    logger << std::setprecision(12) << std::scientific;
    logger << "      " << grid.label() << ".sum = " << sum << "\n";

    for (int z = z0; z < z1; ++z) {
      logger << "        " << z << ", ";

      for (int y = y0; y < y1; ++y) {
        logger << y << ", " << x0 << " | ";

        for (int x = x0; x < x1; ++x) {
          logger << std::setprecision(12) << hgrid(z * gx * gy + y * gx + x)
                 << ", ";
        }

        logger << "\n";
      }
    }
  }

  void
  print_statistics(Bunch& bunch, Logger& logger)
  {

    logger << "Bunch statistics: "
           << "num_valid = " << bunch.get_local_num()
           << ", size = " << bunch.size() << ", capacity = " << bunch.capacity()
           << ", total_num = " << bunch.get_total_num() << "\nMean and std: ";

    // print particles after propagate
    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std = Core_diagnostics::calculate_std(bunch, mean);

    logger << std::setprecision(16) << std::showpos << std::scientific << "\n"
      //<< "\nmean\tstd\n"
      ;

    for (int i = 0; i < 6; ++i) logger << mean[i] << ", " << std[i] << "\n";

    logger << "\n";

    for (int p = 0; p < 4; ++p) bunch.print_particle(p, logger);

    logger << "\n";
  }

  double
  get_smallest_non_tiny(double val, double other1, double other2, double tiny)
  {
    double retval;
    if (val > tiny) {
      retval = val;
    } else {
      if ((other1 > tiny) && (other2 > tiny)) {
        retval = std::min(other1, other2);
      } else {
        retval = std::max(other1, other2);
      }
    }

    return retval;
  }

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

  struct alg_zeroer {
    karray1d_dev arr;

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      arr(i) = 0.0;
    }
  };

  // G000 is naively infinite. In the correct approach, it should be
  // the value which gives the proper integral when convolved with the
  // charge density. Even assuming a constant charge density, the proper
  // value for G000 cannot be computed in closed form. Fortunately,
  // the solver results are insensitive to the exact value of G000.
  // I make the following argument: G000 should be greater than any of
  // the neighboring values of G. The form
  //    G000 = coeff/min(hx,hy,hz),
  // with
  //    coeff > 1
  // satisfies the criterion. An empirical study (see the 3d_open_hockney.py
  // script in docs/devel/solvers) gives coeff = 2.8.
  //
  // const double coeff = 2.8;
  // double G000 = coeff / std::min(hx, std::min(hy, hz));

  // In the following loops we use mirroring for ix and iy, and iz.
  // Note that the doubling algorithm is not quite symmetric. For
  // example, the doubled grid for 4 points in 1d looks like
  //    0 1 2 3 4 3 2 1
  //
  // gx is not the original grid shape in x, instead it is
  //    gx = shape[0] + 1
  // because of the asymmetric in the mirroring

  struct alg_g2_pointlike {
    const double epsilon = 0.01;

    karray1d_dev g2;
    int gx, gy;
    int dgx, dgy, dgz;
    int padded_dgx;
    double hx, hy, hz;

    double igxgy;
    double igx;

    double G000;

    alg_g2_pointlike(karray1d_dev const& g2,
                     std::array<int, 3> const& g,
                     std::array<int, 3> const& dg,
                     std::array<double, 3> const& h)
      : g2(g2)
      , gx(g[0] + 1)
      , gy(g[1] + 1)
      , dgx(dg[0])
      , dgy(dg[1])
      , dgz(dg[2])
      , padded_dgx(Distributed_fft3d::get_padded_shape_real(dgx))
      , hx(h[0])
      , hy(h[1])
      , hz(h[2])
      , igxgy(1.0 / (gx * gy))
      , igx(1.0 / gx)
      , G000(2.8 / std::min(hx, std::min(hy, hz)))
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      int iz = i * igxgy;
      int iy = (i - iz * gx * gy) * igx;
      int ix = i - iz * gx * gy - iy * gx;

      double dx, dy, dz;

      dx = ix * hx;
      dy = iy * hy;
      dz = iz * hz;

      double G;

      if (ix == 0 && iy == 0 && iz == 0)
        G = G000;
      else
        G = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);

#if 0
            if (periodic_z)
            {
                for (int image = -num_images; image <= num_images; ++image)
                {
                    if (image != 0)
                    {
                        double dz_img = dz + image * z_period;
                        const double tiny = 1.0e-9;

                        if (ix==0 && iy==0 && (std::abs(dz_img) < tiny)) G+=G000;
                        else G += 1.0/sqrt(dx*dx + dy*dy + dz_img*dz_img);
                    }
                }
            }
#endif

      int mix, miy, miz;

      mix = dgx - ix;
      if (mix == dgx) mix = 0;

      miy = dgy - iy;
      if (miy == dgy) miy = 0;

      miz = dgz - iz;
      if (miz == dgz) miz = 0;

      g2(iz * padded_dgx * dgy + iy * padded_dgx + ix) = G;
      g2(iz * padded_dgx * dgy + miy * padded_dgx + ix) = G;
      g2(iz * padded_dgx * dgy + iy * padded_dgx + mix) = G;
      g2(iz * padded_dgx * dgy + miy * padded_dgx + mix) = G;

      g2(miz * padded_dgx * dgy + iy * padded_dgx + ix) = G;
      g2(miz * padded_dgx * dgy + miy * padded_dgx + ix) = G;
      g2(miz * padded_dgx * dgy + iy * padded_dgx + mix) = G;
      g2(miz * padded_dgx * dgy + miy * padded_dgx + mix) = G;
    }
  };

  struct alg_g2_linear {
    const double epsilon = 0.01;

    karray1d_dev g2;
    int gx, gy;
    int gx0, gy0, gz0;
    int dgx, dgy, dgz;
    int padded_dgx;
    double hx, hy, hz;

    double igxgy;
    double igx;

    double G000;

    alg_g2_linear(karray1d_dev const& g2,
                  std::array<int, 3> const& g,
                  std::array<int, 3> const& dg,
                  std::array<double, 3> const& h)
      : g2(g2)
      , gx(g[0] + 1)
      , gy(g[1] + 1)
      , gx0(g[0])
      , gy0(g[1])
      , gz0(g[2])
      , dgx(dg[0])
      , dgy(dg[1])
      , dgz(dg[2])
      , padded_dgx(Distributed_fft3d::get_padded_shape_real(dgx))
      , hx(h[0])
      , hy(h[1])
      , hz(h[2])
      , igxgy(1.0 / (gx * gy))
      , igx(1.0 / gx)
      , G000(0.0)
    {
      double rr = hx * hx + hy * hy;
      double r1 = sqrt(hx * hx + hy * hy + hz * hz);
      G000 = (2.0 / rr) * (hz * r1 + rr * log((hz + r1) / sqrt(rr)) -
                           hz * hz); // average value of outer cylinder.
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      int iz = i * igxgy;
      int iy = (i - iz * gx * gy) * igx;
      int ix = i - iz * gx * gy - iy * gx;

      double z = (iz > gz0) ? (dgz - iz) * hz : iz * hz;
      double y = iy * hy;
      double x = ix * hx;

      double rr = x * x + y * y;
      double s_rr_0 = sqrt(rr + z * z);
      double s_rr_n = sqrt(rr + (z - hz) * (z - hz));
      double s_rr_p = sqrt(rr + (z + hz) * (z + hz));

      double G = 2.0 * s_rr_0 - s_rr_n - s_rr_p;

      double T1, T2, r1, r2;
      const double epsz = 1e-12 * hz;

      if (z < -hz) {
        r1 = (s_rr_n - z + hz) / (s_rr_0 - z);
        T1 = (hz - z) * log(r1);
        r2 = (s_rr_0 - z) / (s_rr_p - z - hz);
        T2 = (hz + z) * log(r2);
        G += T1 + T2;
      } else if (std::abs(z + hz) < epsz) {
        r1 = (s_rr_n - z + hz) / (s_rr_0 - z);
        T1 = (hz - z) * log(r1);
        G += T1;
      } else if (std::abs(z) < epsz) {
        if (std::abs(x) + std::abs(y) < 2. * epsz) {
          G += hz * G000;
        } /* T1+T2 in fact */
        else {
          r1 = (sqrt(hz * hz + rr) + hz) / sqrt(rr);
          G += 2.0 * hz * log(r1);
        }
      } else if (std::abs(z - hz) < epsz) {
        r1 = (s_rr_p + z + hz) / (s_rr_0 + z);
        T1 = (hz + z) * log(r1);
        G += T1;
      } else if (z > hz) {
        r1 = (s_rr_0 + z) / (s_rr_n + z - hz);
        T1 = (hz - z) * log(r1);
        r2 = (s_rr_p + z + hz) / (s_rr_0 + z);
        T2 = (hz + z) * log(r2);
        G += T1 + T2;
      } else {
        // throw std::runtime_error(
        //         "Space_charge_3d_open_hockney::get_green_fn2 error1");
      }

      int miy = dgy - iy;
      int mix = dgx - ix;

      if (mix == dgx) mix = 0;
      if (miy == dgy) miy = 0;

      g2(iz * padded_dgx * dgy + iy * padded_dgx + ix) = G;

      if (iy != gy0) {
        g2(iz * padded_dgx * dgy + miy * padded_dgx + ix) = G;

        if (ix != gx0) {
          g2(iz * padded_dgx * dgy + miy * padded_dgx + mix) = G;
        }
      }

      if (ix != gx0) { g2(iz * padded_dgx * dgy + iy * padded_dgx + mix) = G; }
    }
  };

  struct alg_cplx_multiplier {
    karray1d_dev prod;
    karray1d_dev m1, m2;
    int off;

    alg_cplx_multiplier(karray1d_dev const& prod,
                        karray1d_dev const& m1,
                        karray1d_dev const& m2,
                        int offset)
      : prod(prod), m1(m1), m2(m2), off(offset)
    {}

    KOKKOS_INLINE_FUNCTION
    void
    operator()(const int i) const
    {
      const int real = (off + i) * 2;
      const int imag = (off + i) * 2 + 1;

      prod[real] = m1[real] * m2[real] - m1[imag] * m2[imag];
      prod[imag] = m1[real] * m2[imag] + m1[imag] * m2[real];
    }
  };

  struct alg_force_extractor {
    karray1d_dev phi2;
    karray1d_dev enx;
    karray1d_dev eny;
    karray1d_dev enz;

    int gx, gy, gz;
    int dgx, dgy;
    double ihx, ihy, ihz;
    double igxgy;
    double igx;

    alg_force_extractor(karray1d_dev const& phi2,
                        karray1d_dev const& enx,
                        karray1d_dev const& eny,
                        karray1d_dev const& enz,
                        std::array<int, 3> const& g,
                        std::array<int, 3> const& dg,
                        std::array<double, 3> const& h)
      : phi2(phi2)
      , enx(enx)
      , eny(eny)
      , enz(enz)
      , gx(g[0])
      , gy(g[1])
      , gz(g[2])
      , dgx(Distributed_fft3d::get_padded_shape_real(dg[0]))
      , dgy(dg[1])
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

#if 0
            // or, calculate the full domain
            if (ix==0) { ixl = 0; ixr = 1; idx = ihx*2; }
            else if (ix==gx-1) { ixl = gx-2; ixr = gx-1; idx = ihx; }
            else { ixl = ix-1; ixr = ix+1; idx = ihx; }

            if (iy==0) { iyl = 0; iyr = 1; idy = ihy*2; }
            else if (iy==gy-1) { iyl = gy-2; iyr = gy-1; idy = ihy; }
            else { iyl = iy-1; iyr = iy+1; idy = ihy; }

            if (iz==0) { izl = 0; izr = 1; idz = ihz*2; }
            else if (iz==gz-1) { izl = gz-2; izr = gz-1; idz = ihz; }
            else { izl = iz-1; izr = iz+1; idz = ihz; }
#endif

      int idx_r = iz * dgx * dgy + iy * dgx + ixr;
      int idx_l = iz * dgx * dgy + iy * dgx + ixl;
      enx(i) = -(phi2(idx_r) - phi2(idx_l)) * idx;

      idx_r = iz * dgx * dgy + iyr * dgx + ix;
      idx_l = iz * dgx * dgy + iyl * dgx + ix;
      eny(i) = -(phi2(idx_r) - phi2(idx_l)) * idy;

      idx_r = izr * dgx * dgy + iy * dgx + ix;
      idx_l = izl * dgx * dgy + iy * dgx + ix;
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
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
  Space_charge_3d_open_hockney_options const& ops)
  : Collective_operator("sc_3d_open_hockney", 1.0)
  , options(ops)
  , bunch_sim_id()
  , domain(ops.shape, {1.0, 1.0, 1.0})
  , doubled_domain(ops.doubled_shape, {1.0, 1.0, 1.0})
  , use_fixed_domain(false)
  , ffts()
{}

void
Space_charge_3d_open_hockney::apply_impl(Bunch_simulator& sim,
                                         double time_step,
                                         Logger& logger)
{
  logger << "    Space charge 3d open hockney\n";

  scoped_simple_timer timer("sc3d_total");

  // construct the workspace for a new bunch simulator
  if (bunch_sim_id != sim.id()) {
    construct_workspaces(sim);
    bunch_sim_id = sim.id();
  }

  // apply to bunches
  for (size_t t = 0; t < 2; ++t) {
    for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
      apply_bunch(sim[t][b], ffts[t][b], time_step, logger);
    }
  }
}

void
Space_charge_3d_open_hockney::apply_bunch(Bunch& bunch,
                                          Distributed_fft3d& fft,
                                          double time_step,
                                          Logger& logger)
{
  // update domain only when not using fixed
  if (!use_fixed_domain) update_domain(bunch);

  // charge density
  get_local_charge_density(bunch); // [C/m^3]
  get_global_charge_density(bunch);

  // green function
  if (options.green_fn == green_fn_t::pointlike) {
    get_green_fn2_pointlike();
  } else {
    get_green_fn2_linear();
  }

  // potential
  get_local_phi2(fft);
  get_global_phi2(fft);

  auto fn_norm = get_normalization_force(fft);

  // force
  get_force();

  // kick
  apply_kick(bunch, fn_norm, time_step);
}

void
Space_charge_3d_open_hockney::construct_workspaces(Bunch_simulator const& sim)
{
  scoped_simple_timer timer("sc3d_workspaces");

  // doubled shape
  auto const& s = options.doubled_shape;

  // fft objects
  for (size_t t = 0; t < 2; ++t) {
    int num_local_bunches = sim[t].get_bunch_array_size();
    ffts[t] = std::vector<Distributed_fft3d>(num_local_bunches);

    for (size_t b = 0; b < num_local_bunches; ++b) {
      auto comm = sim[t][b].get_comm().divide(options.comm_group_size);

      ffts[t][b].construct(s, comm);
    }
  }

  // local workspaces
  int nx_real = Distributed_fft3d::get_padded_shape_real(s[0]);
  int nx_cplx = Distributed_fft3d::get_padded_shape_cplx(s[0]);

  // doubled domain
  rho2 = karray1d_dev("rho2", nx_real * s[1] * s[2]);
  g2 = karray1d_dev("g2", nx_real * s[1] * s[2]);
  phi2 = karray1d_dev("phi2", nx_real * s[1] * s[2]);

  h_rho2 = Kokkos::create_mirror_view(rho2);
  h_phi2 = Kokkos::create_mirror_view(phi2);

  // En is in the original domain
  enx = karray1d_dev("enx", s[0] * s[1] * s[2] / 8);
  eny = karray1d_dev("eny", s[0] * s[1] * s[2] / 8);
  enz = karray1d_dev("enz", s[0] * s[1] * s[2] / 8);
}

void
Space_charge_3d_open_hockney::update_domain(Bunch const& bunch)
{
  scoped_simple_timer timer("sc3d_domain");

  // do nothing for fixed domain
  if (options.domain_fixed) return;

  auto mean = Core_diagnostics::calculate_mean(bunch);
  auto std = Core_diagnostics::calculate_std(bunch, mean);

  const double tiny = 1.0e-10;

  const auto ix = Bunch::x;
  const auto iy = Bunch::y;
  const auto iz = Bunch::z;

  if ((std[ix] < tiny) && (std[iy] < tiny) && (std[iz] < tiny)) {
    throw std::runtime_error(
      "Space_charge_3d_open_hockney_eigen::update_domain: "
      "all three spatial dimensions have neglible extent");
  }

  std::array<double, 3> offset{mean[ix], mean[iy], mean[iz]};

  std::array<double, 3> size{
    options.n_sigma * get_smallest_non_tiny(std[0], std[2], std[4], tiny),
    options.n_sigma * get_smallest_non_tiny(std[2], std[0], std[4], tiny),
    options.n_sigma * get_smallest_non_tiny(std[4], std[0], std[2], tiny)};

  if (options.grid_entire_period) {
    offset[2] = 0.0;
    size[2] = options.z_period;
  }

  std::array<double, 3> doubled_size{
    size[0] * 2.0, size[1] * 2.0, size[2] * 2.0};

  domain = Rectangular_grid_domain(options.shape, size, offset, false);

  doubled_domain =
    Rectangular_grid_domain(options.doubled_shape, doubled_size, offset, false);
}

void
Space_charge_3d_open_hockney::set_fixed_domain(std::array<double, 3> offset,
                                               std::array<double, 3> size)
{
  if (options.grid_entire_period) {
    offset[2] = 0.0;
    size[2] = options.z_period;
  }

  std::array<double, 3> doubled_size{
    size[0] * 2.0, size[1] * 2.0, size[2] * 2.0};

  domain = Rectangular_grid_domain(options.shape, size, offset, false);

  doubled_domain =
    Rectangular_grid_domain(options.doubled_shape, doubled_size, offset, false);

  use_fixed_domain = true;
}

void
Space_charge_3d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
  scoped_simple_timer timer("sc3d_local_rho");

  auto dg = doubled_domain.get_grid_shape();
  dg[0] = Distributed_fft3d::get_padded_shape_real(dg[0]);

  //#ifdef KOKKOS_ENABLE_CUDA
  deposit_charge_rectangular_3d_kokkos_scatter_view(rho2, domain, dg, bunch);
  //#else
  //  deposit_charge_rectangular_3d_omp_reduce(rho2, domain, dg, bunch);
  //#endif
}

void
Space_charge_3d_open_hockney::get_global_charge_density(Bunch const& bunch)
{
  // do nothing if the bunch occupis a single rank
  if (bunch.get_comm().size() == 1) return;

  scoped_simple_timer timer("sc3d_global_rho");

  auto dg = doubled_domain.get_grid_shape();

  simple_timer_start("sc3d_global_rho_copy");
  Kokkos::deep_copy(h_rho2, rho2);
  simple_timer_stop("sc3d_global_rho_copy");

  simple_timer_start("sc3d_global_rho_reduce");
  int err = MPI_Allreduce(MPI_IN_PLACE,
                          (void*)h_rho2.data(),
                          h_rho2.extent(0),
                          MPI_DOUBLE,
                          MPI_SUM,
                          bunch.get_comm());
  simple_timer_stop("sc3d_global_rho_reduce");

  if (err != MPI_SUCCESS) {
    throw std::runtime_error("MPI error in Space_charge_3d_open_hockney"
                             "(MPI_Allreduce in get_global_charge_density)");
  }

  simple_timer_start("sc3d_global_rho_copy");
  Kokkos::deep_copy(rho2, h_rho2);
  simple_timer_stop("sc3d_global_rho_copy");
}

void
Space_charge_3d_open_hockney::get_green_fn2_pointlike()
{
  if (options.periodic_z) {
    throw std::runtime_error(
      "Space_charge_3d_open_hockney::get_green_fn2_pointlike: "
      "periodic_z not yet implemented");
  }

  scoped_simple_timer timer("sc3d_green_fn2_point");

  auto g = domain.get_grid_shape();
  auto h = doubled_domain.get_cell_size();
  auto dg = doubled_domain.get_grid_shape();

  alg_zeroer az{g2};
  Kokkos::parallel_for(g2.extent(0), az);

  // calculation is performed on grid (gx+1, gy+1, gz+1)
  // rest of the doubled domain will be filled with mirrors
  alg_g2_pointlike alg(g2, g, dg, h);
  Kokkos::parallel_for((g[0] + 1) * (g[1] + 1) * (g[2] + 1), alg);
  Kokkos::fence();
}

void
Space_charge_3d_open_hockney::get_green_fn2_linear()
{
  if (options.periodic_z) {
    throw std::runtime_error(
      "Space_charge_3d_open_hockney::get_green_fn2_linear: "
      "periodic_z not yet implemented");
  }

  scoped_simple_timer timer("sc3d_green_fn2_linear");

  auto g = domain.get_grid_shape();
  auto h = doubled_domain.get_cell_size();
  auto dg = doubled_domain.get_grid_shape();

  alg_zeroer az{g2};
  Kokkos::parallel_for(g2.extent(0), az);

  // calculation is performed on grid (gx+1, gy+1, gz+1)
  // rest of the doubled domain will be filled with mirrors
  alg_g2_linear alg(g2, g, dg, h);
  Kokkos::parallel_for((g[0] + 1) * (g[1] + 1) * dg[2], alg);
  Kokkos::fence();
}

void
Space_charge_3d_open_hockney::get_local_phi2(Distributed_fft3d& fft)
{
  scoped_simple_timer timer("sc3d_local_f");

  auto dg = doubled_domain.get_grid_shape();

  // FFT
  fft.transform(rho2, rho2);
  Kokkos::fence();

  fft.transform(g2, g2);
  Kokkos::fence();

  // zero phi2 when using multiple ranks
  if (fft.get_comm().size() > 1) {
    int padded_gx_real = fft.padded_nx_real();
    alg_zeroer az{phi2};
    Kokkos::parallel_for(padded_gx_real * dg[0] * dg[1], az);
  }

  int lower = fft.get_lower();
  int upper = fft.get_upper();

  int padded_gx_cplx = fft.padded_nx_cplx();
  int offset = lower * padded_gx_cplx * dg[1];
  int nz = upper - lower;

  alg_cplx_multiplier alg(phi2, rho2, g2, offset);
  Kokkos::parallel_for(nz * dg[1] * padded_gx_cplx, alg);
  Kokkos::fence();

  // inv fft
  fft.inv_transform(phi2, phi2);
  Kokkos::fence();
}

void
Space_charge_3d_open_hockney::get_global_phi2(Distributed_fft3d const& fft)
{
  // do nothing if the solver only has a single rank
  if (fft.get_comm().size() == 1) return;

  scoped_simple_timer timer("sc3d_global_f");

  Kokkos::deep_copy(h_phi2, phi2);

  auto dg = doubled_domain.get_grid_shape();
  auto nx_real = fft.padded_nx_real();

  int err = MPI_Allreduce(MPI_IN_PLACE,
                          (void*)h_phi2.data(),
                          nx_real * dg[1] * dg[2],
                          MPI_DOUBLE,
                          MPI_SUM,
                          fft.get_comm());

  if (err != MPI_SUCCESS) {
    throw std::runtime_error(
      "MPI error in Space_charge_3d_open_hockney"
      "(MPI_Allreduce in get_global_electric_force2_allreduce)");
  }

  Kokkos::deep_copy(phi2, h_phi2);
}

double
Space_charge_3d_open_hockney::get_normalization_force(
  Distributed_fft3d const& fft)
{
  auto h = domain.get_cell_size();

  double hx = h[0];
  double hy = h[1];
  double hz = h[2];

  // volume element in integral
  double normalization = hx * hy * hz;

  // dummy factor from weight0 of deposit.cc
  normalization *= 1.0 / (4.0 * pi * pconstants::epsilon0);

  // from charege density
  normalization *= 1.0;

  // 1.0 from point-like greens function.
  // 1.0/(hz*hz) for linear greens function
  if (options.green_fn == green_fn_t::linear) normalization *= 1.0 / (hz * hz);

  // from fft
  normalization *= fft.get_roundtrip_normalization();

  return normalization;
}

void
Space_charge_3d_open_hockney::get_force()
{
  auto g = domain.get_grid_shape();
  auto h = doubled_domain.get_cell_size();
  auto dg = doubled_domain.get_grid_shape();

  // phi2 is in (padded_real_dgx, dgy, dgz)
  // en{x|y|z} is in (gx, gy, gz)
  alg_force_extractor alg(phi2, enx, eny, enz, g, dg, h);
  Kokkos::parallel_for(g[0] * g[1] * g[2], alg);
  Kokkos::fence();
}

void
Space_charge_3d_open_hockney::apply_kick(Bunch& bunch,
                                         double fn_norm,
                                         double time_step)
{
  scoped_simple_timer timer("sc3d_kick");

  auto ref = bunch.get_reference_particle();

  double q = bunch.get_particle_charge() * pconstants::e;
  double m = bunch.get_mass();

  double gamma = ref.get_gamma();
  double beta = ref.get_beta();
  double pref = ref.get_momentum();

  double unit_conversion = pconstants::c / (1e9 * pconstants::e);
  double factor = options.kick_scale * unit_conversion * q * time_step *
                  fn_norm / (pref * gamma * gamma * beta);

  auto parts = bunch.get_local_particles();
  auto masks = bunch.get_local_particle_masks();

  auto g = domain.get_grid_shape();
  auto h = domain.get_cell_size();
  auto l = domain.get_left();

  alg_kicker kicker(parts, masks, enx, eny, enz, g, h, l, factor, pref, m);

  Kokkos::parallel_for(bunch.size(), kicker);
  Kokkos::fence();
}
