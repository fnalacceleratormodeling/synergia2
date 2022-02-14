
#include "space_charge_3d_open_hockney.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"
#include "deposit.h"

#include "synergia/utils/simple_timer.h"

using mconstants::pi;
using pconstants::epsilon0;

namespace
{
    void
    print_grid( Logger & logger, 
                karray1d_dev const & grid, 
                int x0, int x1, 
                int y0, int y1,
                int z0, int z1,
                int gx, int gy, int gz,
                int off = 0 )
    {
        karray1d hgrid = Kokkos::create_mirror_view(grid);
        Kokkos::deep_copy(hgrid, grid);

        double sum = 0;

        int dim = grid.extent(0);
        for(int i=0; i<dim; ++i) 
        {
            sum += fabs(hgrid(i));
        }

#if 0
        for(int x=0; x<gx; ++x)
            for(int y=0; y<gy; ++y)
                sum += hgrid((x*gy + y)*2 + off);
#endif
        
        logger << std::setprecision(12) << std::scientific;
        logger << "      " << grid.label() << ".sum = " << sum << "\n";

        for (int z=z0; z<z1; ++z)
        {
            logger << "        " << z << ", ";

            for (int y=y0; y<y1; ++y)
            {
                logger << y << ", " << x0 << " | ";

                for (int x=x0; x<x1; ++x)
                {
                    logger << std::setprecision(12) 
                        << hgrid(z*gx*gy + y*gx + x) << ", ";
                }

                logger << "\n";
            }
        }
    }

    void print_statistics(Bunch & bunch, Logger & logger)
    {

        logger
            << "Bunch statistics: "
            << "num_valid = " << bunch.get_local_num()
            << ", size = " << bunch.size()
            << ", capacity = " << bunch.capacity()
            << ", total_num = " << bunch.get_total_num()
            <<"\nMean and std: ";


        // print particles after propagate
        auto mean = Core_diagnostics::calculate_mean(bunch);
        auto std  = Core_diagnostics::calculate_std(bunch, mean);

        logger
            << std::setprecision(16)
            << std::showpos << std::scientific
            << "\n"
            //<< "\nmean\tstd\n"
            ;

        for (int i=0; i<6; ++i) 
            logger << mean[i] << ", " << std[i] << "\n";

        logger << "\n";

        for (int p=0; p<4; ++p) bunch.print_particle(p, logger);

        logger << "\n";
    }

    double
    get_smallest_non_tiny( double val, 
                           double other1, 
                           double other2, 
                           double tiny )
    {
        double retval;
        if (val > tiny) 
        {
            retval = val;
        } 
        else 
        {
            if ((other1 > tiny) && (other2 > tiny)) 
            {
                retval = std::min(other1, other2);
            } 
            else 
            {
                retval = std::max(other1, other2);
            }
        }

        return retval;
    }

    KOKKOS_INLINE_FUNCTION
    int fast_int_floor_kokkos(const double x)
    {
        int ix = static_cast<int>(x);
        return x > 0.0 ? ix : ((x - ix == 0) ? ix : ix - 1);
    }

    KOKKOS_INLINE_FUNCTION
    void get_leftmost_indices_offset(
            double pos, double left, double inv_cell_size,
            int & idx, double & off )
    {
        double scaled_location = (pos - left) * inv_cell_size - 0.5;
        idx = fast_int_floor_kokkos(scaled_location);
        off = scaled_location - idx;
    }

    struct alg_zeroer
    {
        karray1d_dev arr;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { arr(i) = 0.0; }
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

    struct alg_g2_pointlike
    {
        const double epsilon = 0.01;

        karray1d_dev g2;
        int gx, gy;
        int dgx, dgy, dgz;
        int padded_dgx;
        double hx, hy, hz;

        double igxgy;
        double igx;

        double G000;

        alg_g2_pointlike(
                karray1d_dev const & g2,
                std::array<int, 3> const & g,
                std::array<int, 3> const & dg,
                std::array<double, 3> const & h ) 
            : g2(g2)
            , gx(g[0]+1), gy(g[1]+1)
            , dgx(dg[0]), dgy(dg[1]), dgz(dg[2])
            , padded_dgx(Distributed_fft3d::get_padded_shape_real(dgx))
            , hx(h[0]), hy(h[1]), hz(h[2])
            , igxgy(1.0/(gx*gy))
            , igx(1.0/gx)
            , G000(2.8/std::min(hx, std::min(hy, hz)))
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz*gx*gy) * igx;
            int ix = i - iz*gx*gy - iy*gx;

            double dx, dy, dz;

            dx = ix * hx;
            dy = iy * hy;
            dz = iz * hz;

            double G;

            if (ix==0 && iy==0 && iz==0) G = G000;
            else G = 1.0 / sqrt(dx*dx + dy*dy + dz*dz);

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

            g2(iz*padded_dgx*dgy +  iy*padded_dgx +  ix) = G;
            g2(iz*padded_dgx*dgy + miy*padded_dgx +  ix) = G;
            g2(iz*padded_dgx*dgy +  iy*padded_dgx + mix) = G;
            g2(iz*padded_dgx*dgy + miy*padded_dgx + mix) = G;

            g2(miz*padded_dgx*dgy +  iy*padded_dgx +  ix) = G;
            g2(miz*padded_dgx*dgy + miy*padded_dgx +  ix) = G;
            g2(miz*padded_dgx*dgy +  iy*padded_dgx + mix) = G;
            g2(miz*padded_dgx*dgy + miy*padded_dgx + mix) = G;
        }
    };

    struct alg_g2_linear
    {
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

        alg_g2_linear(
                karray1d_dev const & g2,
                std::array<int, 3> const & g,
                std::array<int, 3> const & dg,
                std::array<double, 3> const & h ) 
            : g2(g2)
            , gx(g[0]+1), gy(g[1]+1)
            , gx0(g[0]), gy0(g[1]), gz0(g[2])
            , dgx(dg[0]), dgy(dg[1]), dgz(dg[2])
            , padded_dgx(Distributed_fft3d::get_padded_shape_real(dgx))
            , hx(h[0]), hy(h[1]), hz(h[2])
            , igxgy(1.0/(gx*gy))
            , igx(1.0/gx)
            , G000(0.0)
        { 
            double rr = hx * hx + hy * hy;
            double r1 = sqrt(hx * hx + hy * hy + hz * hz);
            G000 = (2.0 / rr) * (hz * r1 + rr * log((hz + r1) 
                        / sqrt(rr)) - hz * hz); // average value of outer cylinder.
        }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz*gx*gy) * igx;
            int ix = i - iz*gx*gy - iy*gx;

            double z = (iz>gz0) ? (dgz-iz)*hz : iz*hz;
            double y = iy * hy;
            double x = ix * hx;

            double rr = x * x + y * y;
            double s_rr_0 = sqrt(rr + z*z);
            double s_rr_n = sqrt(rr + (z-hz)*(z-hz));
            double s_rr_p = sqrt(rr + (z+hz)*(z+hz));

            double G = 2.0 * s_rr_0 - s_rr_n - s_rr_p;

            double T1, T2, r1, r2;
            const double epsz = 1e-12 * hz;

            if (z < -hz) 
            {
                r1 = (s_rr_n - z + hz) / (s_rr_0 - z);
                T1 = (hz - z) * log(r1);
                r2 = (s_rr_0 - z) / (s_rr_p - z - hz);
                T2 = (hz + z) * log(r2);
                G += T1 + T2;
            } 
            else if (std::abs(z + hz) < epsz) 
            {
                r1 = (s_rr_n - z + hz) / (s_rr_0 - z);
                T1 = (hz - z) * log(r1);
                G += T1;
            } 
            else if (std::abs(z) < epsz) 
            {
                if (std::abs(x) + std::abs(y) < 2. * epsz) 
                {
                    G += hz * G000;
                } /* T1+T2 in fact */
                else 
                {
                    r1 = (sqrt(hz * hz + rr) + hz) / sqrt(rr);
                    G += 2.0 * hz * log(r1);
                }
            } 
            else if (std::abs(z - hz) < epsz) 
            {
                r1 = (s_rr_p + z + hz) / (s_rr_0 + z);
                T1 = (hz + z) * log(r1);
                G += T1;
            } 
            else if (z > hz) 
            {
                r1 = (s_rr_0 + z) / (s_rr_n + z - hz);
                T1 = (hz - z) * log(r1);
                r2 = (s_rr_p + z + hz) / (s_rr_0 + z);
                T2 = (hz + z) * log(r2);
                G += T1 + T2;
            } 
            else 
            {
                //throw std::runtime_error(
                //        "Space_charge_3d_open_hockney::get_green_fn2 error1");
            }

            int miy = dgy - iy;
            int mix = dgx - ix;

            if (mix == dgx) mix = 0;
            if (miy == dgy) miy = 0;

            g2(iz*padded_dgx*dgy + iy*padded_dgx + ix) = G;

            if (iy != gy0)
            {
                g2(iz*padded_dgx*dgy + miy*padded_dgx + ix) = G;

                if (ix != gx0)
                {
                    g2(iz*padded_dgx*dgy + miy*padded_dgx + mix) = G;
                }
            }

            if (ix != gx0)
            {
                g2(iz*padded_dgx*dgy + iy*padded_dgx + mix) = G;
            }
        }
    };

    struct alg_cplx_multiplier
    {
        karray1d_dev prod;
        karray1d_dev m1, m2;
        int off;

        alg_cplx_multiplier(
                karray1d_dev const & prod,
                karray1d_dev const & m1,
                karray1d_dev const & m2,
                int offset )
            : prod(prod), m1(m1), m2(m2), off(offset)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            const int real = (off+i)*2;
            const int imag = (off+i)*2 + 1;

            prod[real] = m1[real]*m2[real] - m1[imag]*m2[imag];
            prod[imag] = m1[real]*m2[imag] + m1[imag]*m2[real];
        }
    };

    struct alg_force_extractor
    {
        karray1d_dev phi2;
        karray1d_dev enx;
        karray1d_dev eny;
        karray1d_dev enz;

        int gx, gy, gz;
        int dgx, dgy;
        double ihx, ihy, ihz;
        double igxgy;
        double igx;

        alg_force_extractor(
                karray1d_dev const& phi2,
                karray1d_dev const& enx,
                karray1d_dev const& eny,
                karray1d_dev const& enz,
                std::array<int, 3> const& g,
                std::array<int, 3> const& dg,
                std::array<double, 3> const& h )
            : phi2(phi2), enx(enx), eny(eny), enz(enz)
            , gx(g[0]), gy(g[1]), gz(g[2])
            , dgx(Distributed_fft3d::get_padded_shape_real(dg[0]))
            , dgy(dg[1])
            , ihx(0.5/h[0]), ihy(0.5/h[1]), ihz(0.5/h[2])
            , igxgy(1.0/(gx*gy))
            , igx(1.0/gx)
        { }


        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int iz = i * igxgy;
            int iy = (i - iz*gx*gy) * igx;
            int ix = i - iz*gx*gy - iy*gx;

            int ixl, ixr, iyl, iyr, izl, izr;

            double idx, idy, idz;

            // all boundaries will be skipped (set to 0)
            if ( ix==0 || ix==gx-1 || 
                 iy==0 || iy==gy-1 || 
                 iz==0 || iz==gz-1    )
                return;

            ixl = ix - 1; ixr = ix + 1;
            iyl = iy - 1; iyr = iy + 1;
            izl = iz - 1; izr = iz + 1;

            idx = ihx; idy = ihy; idz = ihz;

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

            int idx_r = iz*dgx*dgy + iy*dgx + ixr;
            int idx_l = iz*dgx*dgy + iy*dgx + ixl;
            enx(i) = -(phi2(idx_r) - phi2(idx_l)) * idx;

            idx_r = iz*dgx*dgy + iyr*dgx + ix;
            idx_l = iz*dgx*dgy + iyl*dgx + ix;
            eny(i) = -(phi2(idx_r) - phi2(idx_l)) * idy;

            idx_r = izr*dgx*dgy + iy*dgx + ix;
            idx_l = izl*dgx*dgy + iy*dgx + ix;
            enz(i) = -(phi2(idx_r) - phi2(idx_l)) * idz;
        }
    };

    struct alg_kicker
    {
        Particles parts;
        ConstParticleMasks masks;

        karray1d_dev enx;
        karray1d_dev eny;
        karray1d_dev enz;

        int gx, gy, gz;
        double ihx, ihy, ihz;
        double lx, ly, lz;
        double factor, pref, m;

        alg_kicker(
                Particles parts,
                ConstParticleMasks masks,
                karray1d_dev const& enx,
                karray1d_dev const& eny,
                karray1d_dev const& enz,
                std::array<int,    3> const& g,
                std::array<double, 3> const& h,
                std::array<double, 3> const& l,
                double factor,
                double pref,
                double m )
            : parts(parts), masks(masks)
            , enx(enx), eny(eny), enz(enz)
            , gx(g[0]), gy(g[1]), gz(g[2])
            , ihx(1.0/h[0]), ihy(1.0/h[1]), ihz(1.0/h[2])
            , lx(l[0]), ly(l[1]), lz(l[2])
            , factor(factor), pref(pref), m(m)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            if (masks(i))
            {
                int ix, iy, iz;
                double ox, oy, oz;

                get_leftmost_indices_offset(parts(i, 0), lx, ihx, ix, ox);
                get_leftmost_indices_offset(parts(i, 2), ly, ihy, iy, oy);
                get_leftmost_indices_offset(parts(i, 4), lz, ihz, iz, oz);

                double aox = 1.0 - ox;
                double aoy = 1.0 - oy;
                double aoz = 1.0 - oz;

                if ( (ix>=0 && ix<gx-1) && 
                     (iy>=0 && iy<gy-1) && 
                     ((iz>=0 && iz<gz-1) /* || periodic_z */) ) 
                {
                    double val = 0;
                    int base = 0;

                    // enz
                    base = iz*gx*gy + iy*gx + ix;
                    val  = aoz * aoy * aox * enz(base);      // z, y, x
                    val += aoz * aoy *  ox * enz(base+1);    // z, y, x+1
                    val += aoz *  oy * aox * enz(base+gx);   // z, y+1, x
                    val += aoz *  oy *  ox * enz(base+gx+1); // z, y+1, x+1

                    base = (iz+1)*gx*gy + iy*gx + ix;
                    val += oz * aoy * aox * enz(base);       // z+1, y, x
                    val += oz * aoy *  ox * enz(base+1);     // z+1, y, x+1
                    val += oz *  oy * aox * enz(base+gx);    // z+1, y+1, x
                    val += oz *  oy *  ox * enz(base+gx+1);  // z+1, y+1, x+1

                    double p = pref + parts(i, 5) * pref;
                    double Eoc_i = std::sqrt(p*p+m*m);
                    double Eoc_f = Eoc_i + factor * (-pref) * val;
                    double dpop = (std::sqrt(Eoc_f*Eoc_f-m*m) - std::sqrt(Eoc_i*Eoc_i-m*m))/pref;

                    parts(i, 5) += dpop;

                    // eny
                    base = iz*gx*gy + iy*gx + ix;
                    val  = aoz * aoy * aox * eny(base);      // z, y, x
                    val += aoz * aoy *  ox * eny(base+1);    // z, y, x+1
                    val += aoz *  oy * aox * eny(base+gx);   // z, y+1, x
                    val += aoz *  oy *  ox * eny(base+gx+1); // z, y+1, x+1

                    base = (iz+1)*gx*gy + iy*gx + ix;
                    val += oz * aoy * aox * eny(base);       // z+1, y, x
                    val += oz * aoy *  ox * eny(base+1);     // z+1, y, x+1
                    val += oz *  oy * aox * eny(base+gx);    // z+1, y+1, x
                    val += oz *  oy *  ox * eny(base+gx+1);  // z+1, y+1, x+1

                    parts(i, 3) += factor * val;

                    // enx
                    base = iz*gx*gy + iy*gx + ix;
                    val  = aoz * aoy * aox * enx(base);      // z, y, x
                    val += aoz * aoy *  ox * enx(base+1);    // z, y, x+1
                    val += aoz *  oy * aox * enx(base+gx);   // z, y+1, x
                    val += aoz *  oy *  ox * enx(base+gx+1); // z, y+1, x+1

                    base = (iz+1)*gx*gy + iy*gx + ix;
                    val += oz * aoy * aox * enx(base);       // z+1, y, x
                    val += oz * aoy *  ox * enx(base+1);     // z+1, y, x+1
                    val += oz *  oy * aox * enx(base+gx);    // z+1, y+1, x
                    val += oz *  oy *  ox * enx(base+gx+1);  // z+1, y+1, x+1

                    parts(i, 1) += factor * val;
                }
            }
        }
    };
}


Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        Space_charge_3d_open_hockney_options const & ops)
    : Collective_operator("sc_3d_open_hockney", 1.0)
    , options(ops)
    , bunch_sim_id()
    , domain(ops.shape, {1.0, 1.0, 1.0})
    , doubled_domain(ops.doubled_shape, {1.0, 1.0, 1.0})
    , use_fixed_domain(false)
    , ffts()
{
}



void 
Space_charge_3d_open_hockney::apply_impl(
            Bunch_simulator& sim, 
            double time_step, 
            Logger& logger)
{
    logger << "    Space charge 3d open hockney\n";

    scoped_simple_timer timer("sc3d_total");

    // construct the workspace for a new bunch simulator
    if (bunch_sim_id != sim.id())
    {
        construct_workspaces(sim);
        bunch_sim_id = sim.id();
    }

    // apply to bunches
    for(size_t t=0; t<2; ++t)
    {
        for(size_t b=0; b<sim[t].get_bunch_array_size(); ++b)
        {
            apply_bunch(sim[t][b], ffts[t][b], time_step, logger);
        }
    }
}

void
Space_charge_3d_open_hockney::apply_bunch(
            Bunch& bunch, 
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
    if (options.green_fn == green_fn_t::pointlike) 
    {
        get_green_fn2_pointlike();
    }
    else 
    {
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
Space_charge_3d_open_hockney::construct_workspaces(
        Bunch_simulator const& sim)
{
    scoped_simple_timer timer("sc3d_workspaces");

    // doubled shape
    auto const& s = options.doubled_shape;

    // fft objects
    for(size_t t=0; t<2; ++t)
    {
        int num_local_bunches = sim[t].get_bunch_array_size();
        ffts[t] = std::vector<Distributed_fft3d>(num_local_bunches);

        for(size_t b=0; b<num_local_bunches; ++b)
        {
            auto comm = sim[t][b]
                .get_comm()
                .divide(options.comm_group_size);

            ffts[t][b].construct(s, comm);
        }
    }

    // local workspaces
    int nx_real = Distributed_fft3d::get_padded_shape_real(s[0]);
    int nx_cplx = Distributed_fft3d::get_padded_shape_cplx(s[0]);

    // doubled domain
    rho2 = karray1d_dev("rho2", nx_real * s[1] * s[2]);
      g2 = karray1d_dev(  "g2", nx_real * s[1] * s[2]);
    phi2 = karray1d_dev("phi2", nx_real * s[1] * s[2]);

    h_rho2 = Kokkos::create_mirror_view(rho2);
    h_phi2 = Kokkos::create_mirror_view(phi2);

    // En is in the original domain
    enx = karray1d_dev("enx", s[0] * s[1] * s[2] / 8);
    eny = karray1d_dev("eny", s[0] * s[1] * s[2] / 8);
    enz = karray1d_dev("enz", s[0] * s[1] * s[2] / 8);
}

void
Space_charge_3d_open_hockney::update_domain(Bunch const & bunch)
{
    scoped_simple_timer timer("sc3d_domain");

    // do nothing for fixed domain
    if (options.domain_fixed) return;

    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std  = Core_diagnostics::calculate_std(bunch, mean);

    const double tiny = 1.0e-10;

    const auto ix = Bunch::x;
    const auto iy = Bunch::y;
    const auto iz = Bunch::z;

    if ( (std[ix] < tiny) && (std[iy] < tiny) && (std[iz] < tiny) ) 
    {
        throw std::runtime_error(
                "Space_charge_3d_open_hockney_eigen::update_domain: "
                "all three spatial dimensions have neglible extent");
    }

    std::array<double, 3> offset { 
        mean[ix], 
        mean[iy], 
        mean[iz] };

    std::array<double, 3> size { 
        options.n_sigma * get_smallest_non_tiny(std[0], std[2], std[4], tiny),
        options.n_sigma * get_smallest_non_tiny(std[2], std[0], std[4], tiny),
        options.n_sigma * get_smallest_non_tiny(std[4], std[0], std[2], tiny) };

    if (options.grid_entire_period)
    {
        offset[2] = 0.0;
        size[2] = options.z_period;
    }

    std::array<double, 3> doubled_size { 
        size[0] * 2.0, 
        size[1] * 2.0, 
        size[2] * 2.0 };

    domain = Rectangular_grid_domain(
            options.shape, size, offset, false);

    doubled_domain = Rectangular_grid_domain(
            options.doubled_shape, doubled_size, offset, false);
}

void
Space_charge_3d_open_hockney::set_fixed_domain(
        std::array<double, 3> offset, 
        std::array<double, 3> size)
{
    if (options.grid_entire_period)
    {
        offset[2] = 0.0;
        size[2] = options.z_period;
    }

    std::array<double, 3> doubled_size { 
        size[0] * 2.0, 
        size[1] * 2.0, 
        size[2] * 2.0 };

    domain = Rectangular_grid_domain(
            options.shape, size, offset, false);

    doubled_domain = Rectangular_grid_domain(
            options.doubled_shape, doubled_size, offset, false);

    use_fixed_domain = true;
}


void
Space_charge_3d_open_hockney::get_local_charge_density(
        Bunch const& bunch)
{
    scoped_simple_timer timer("sc3d_local_rho");

    auto dg = doubled_domain.get_grid_shape();
    dg[0] = Distributed_fft3d::get_padded_shape_real(dg[0]);

#ifdef KOKKOS_ENABLE_CUDA
    deposit_charge_rectangular_3d_kokkos_scatter_view(rho2,
            domain, dg, bunch);
#else
    deposit_charge_rectangular_3d_omp_reduce(rho2,
            domain, dg, bunch);
#endif
}

void
Space_charge_3d_open_hockney::get_global_charge_density(
        Bunch const & bunch )
{
    // do nothing if the bunch occupis a single rank
    if (bunch.get_comm().size() == 1) return;

    scoped_simple_timer timer("sc3d_global_rho");

    auto dg = doubled_domain.get_grid_shape();

    simple_timer_start("sc3d_global_rho_copy");
    Kokkos::deep_copy(h_rho2, rho2);
    simple_timer_stop("sc3d_global_rho_copy");

    simple_timer_start("sc3d_global_rho_reduce");
    int err = MPI_Allreduce( MPI_IN_PLACE,
                             (void*)h_rho2.data(), 
                             h_rho2.extent(0),
                             MPI_DOUBLE, 
                             MPI_SUM, 
                             bunch.get_comm() );
    simple_timer_stop("sc3d_global_rho_reduce");

    if (err != MPI_SUCCESS)
    {
        throw std::runtime_error( 
                "MPI error in Space_charge_3d_open_hockney"
                "(MPI_Allreduce in get_global_charge_density)" );
    }

    simple_timer_start("sc3d_global_rho_copy");
    Kokkos::deep_copy(rho2, h_rho2);
    simple_timer_stop("sc3d_global_rho_copy");
}

void
Space_charge_3d_open_hockney::get_green_fn2_pointlike()
{
    if (options.periodic_z)
    {
        throw std::runtime_error(
                "Space_charge_3d_open_hockney::get_green_fn2_pointlike: "
                "periodic_z not yet implemented");
    }

    scoped_simple_timer timer("sc3d_green_fn2_point");

    auto  g = domain.get_grid_shape();
    auto  h = doubled_domain.get_cell_size();
    auto dg = doubled_domain.get_grid_shape();

    alg_zeroer az{g2};
    Kokkos::parallel_for(g2.extent(0), az);

    // calculation is performed on grid (gx+1, gy+1, gz+1)
    // rest of the doubled domain will be filled with mirrors
    alg_g2_pointlike alg(g2, g, dg, h);
    Kokkos::parallel_for((g[0]+1)*(g[1]+1)*(g[2]+1), alg);
    Kokkos::fence();
}

void
Space_charge_3d_open_hockney::get_green_fn2_linear()
{
    if (options.periodic_z)
    {
        throw std::runtime_error(
                "Space_charge_3d_open_hockney::get_green_fn2_linear: "
                "periodic_z not yet implemented");
    }

    scoped_simple_timer timer("sc3d_green_fn2_linear");

    auto  g = domain.get_grid_shape();
    auto  h = doubled_domain.get_cell_size();
    auto dg = doubled_domain.get_grid_shape();

    alg_zeroer az{g2};
    Kokkos::parallel_for(g2.extent(0), az);

    // calculation is performed on grid (gx+1, gy+1, gz+1)
    // rest of the doubled domain will be filled with mirrors
    alg_g2_linear alg(g2, g, dg, h);
    Kokkos::parallel_for((g[0]+1)*(g[1]+1)*dg[2], alg);
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

    fft.transform(  g2,   g2);
    Kokkos::fence();

    // zero phi2 when using multiple ranks
    if (fft.get_comm().size() > 1)
    {
        int padded_gx_real = fft.padded_nx_real();
        alg_zeroer az{phi2};
        Kokkos::parallel_for(padded_gx_real*dg[0]*dg[1], az);
    }

    int lower = fft.get_lower();
    int upper = fft.get_upper();

    int padded_gx_cplx = fft.padded_nx_cplx();
    int offset = lower*padded_gx_cplx*dg[1];
    int nz = upper - lower;

    alg_cplx_multiplier alg(phi2, rho2, g2, offset);
    Kokkos::parallel_for(nz*dg[1]*padded_gx_cplx, alg);
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

    int err = MPI_Allreduce( MPI_IN_PLACE,
                             (void*)h_phi2.data(), 
                             nx_real*dg[1]*dg[2], 
                             MPI_DOUBLE, 
                             MPI_SUM, 
                             fft.get_comm() );

    if (err != MPI_SUCCESS)
    {
        throw std::runtime_error( 
                "MPI error in Space_charge_3d_open_hockney"
                "(MPI_Allreduce in get_global_electric_force2_allreduce)" );
    }

    Kokkos::deep_copy(phi2, h_phi2);
}

double
Space_charge_3d_open_hockney::get_normalization_force(Distributed_fft3d const& fft)
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
    if (options.green_fn == green_fn_t::linear)
        normalization *= 1.0 / (hz*hz);

    // from fft
    normalization *= fft.get_roundtrip_normalization();

    return normalization;
}

void
Space_charge_3d_open_hockney::get_force()
{
    auto  g = domain.get_grid_shape();
    auto  h = doubled_domain.get_cell_size();
    auto dg = doubled_domain.get_grid_shape();

    // phi2 is in (padded_real_dgx, dgy, dgz)
    // en{x|y|z} is in (gx, gy, gz)
    alg_force_extractor alg(phi2, enx, eny, enz, g, dg, h);
    Kokkos::parallel_for(g[0]*g[1]*g[2], alg);
    Kokkos::fence();
}


void
Space_charge_3d_open_hockney::apply_kick(
        Bunch & bunch,
        double fn_norm,
        double time_step )
{
    scoped_simple_timer timer("sc3d_kick");

    auto ref = bunch.get_reference_particle();

    double q = bunch.get_particle_charge() * pconstants::e;
    double m = bunch.get_mass();

    double gamma = ref.get_gamma();
    double beta  = ref.get_beta();
    double pref  = ref.get_momentum();

    double unit_conversion = pconstants::c / (1e9 * pconstants::e);
    double factor = options.kick_scale * unit_conversion * q * time_step * fn_norm 
        / (pref * gamma * gamma * beta);

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    alg_kicker kicker(parts, masks, enx, eny, enz,
            g, h, l, factor, pref, m);

    Kokkos::parallel_for(bunch.size(), kicker);
    Kokkos::fence();
}





 


#if 0
void
Space_charge_3d_open_hockney::setup_communication(
        Commxx_sptr const& bunch_comm_sptr)
{
    if (comm2_sptr != commxx_divider_sptr->get_commxx_sptr(bunch_comm_sptr)) {
        comm2_sptr = commxx_divider_sptr->get_commxx_sptr(bunch_comm_sptr);
        setup_derived_communication();
    }
}

void
Space_charge_3d_open_hockney::setup_derived_communication()
{
    distributed_fft3d_sptr = Distributed_fft3d_sptr(
            new Distributed_fft3d(doubled_grid_shape, comm2_sptr));
    padded_grid_shape = distributed_fft3d_sptr->get_padded_shape_real();
    std::vector<int > ranks1; // ranks with data from the undoubled domain
    int lower = 0;
    for (int rank = 0; rank < comm2_sptr->get_size(); ++rank) {
        int uppers2 = distributed_fft3d_sptr->get_uppers()[rank];
        int uppers1 = std::min(uppers2, grid_shape[0]);
        int length0;
        if (rank > 0) {
            length0 = uppers1 - distributed_fft3d_sptr->get_uppers()[rank - 1];
        } else {
            length0 = uppers1;
        }
        if (length0 > 0) {
            ranks1.push_back(rank);
            lowers1.push_back(lower);
            int total_length = length0 * grid_shape[1] * grid_shape[2];
            lengths1.push_back(total_length);
            lower += total_length;
        }
    }
    comm1_sptr = Commxx_sptr(new Commxx(comm2_sptr, ranks1));
    std::vector<int > real_uppers(distributed_fft3d_sptr->get_uppers());
    real_lengths = distributed_fft3d_sptr->get_lengths();
    for (int i = 0; i < comm2_sptr->get_size(); ++i) {
        if (real_uppers[i] > grid_shape[0]) {
            real_uppers[i] = grid_shape[0];
        }
        if (i == 0) {
            real_lengths[0] = real_uppers[0] * grid_shape[1] * grid_shape[2];
        } else {
            real_lengths[i] = (real_uppers[i] - real_uppers[i - 1])
                    * grid_shape[1] * grid_shape[2];
        }
    }
    int my_rank = comm2_sptr->get_rank();
    if (my_rank > 0) {
        real_lower = real_uppers[my_rank - 1];
    } else {
        real_lower = 0;
    }
    real_upper = real_uppers[my_rank];
    real_length = real_lengths[my_rank];
    if (my_rank > 0) {
        doubled_lower = distributed_fft3d_sptr->get_uppers()[my_rank - 1];
    } else {
        doubled_lower = 0;
    }
    doubled_upper = distributed_fft3d_sptr->get_uppers()[my_rank];
    real_doubled_lower = std::min(doubled_lower, grid_shape[0]);
    real_doubled_upper = std::min(doubled_upper, grid_shape[0]);
}

void
Space_charge_3d_open_hockney::constructor_common(
        std::vector<int > const& grid_shape)
{
    if (this->periodic_z && (this->z_period == 0.0)) {
        throw std::runtime_error(
                "Space_charge_3d_open_hockney: z_period cannot be 0 when periodic_z is true");
    }
    this->grid_shape[0] = grid_shape[2];
    this->grid_shape[1] = grid_shape[1];
    this->grid_shape[2] = grid_shape[0];
    for (int i = 0; i < 3; ++i) {
        doubled_grid_shape[i] = 2 * this->grid_shape[i];
    }

    set_green_fn_type(linear);
    set_charge_density_comm(charge_allreduce);
    set_e_field_comm(e_field_allreduce);
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        std::vector<int > const & grid_shape, bool longitudinal_kicks,
        bool periodic_z, double z_period, bool grid_entire_period,
        double n_sigma, double kick_scale) :
                Collective_operator("space charge 3D open hockney"),
                grid_shape(3),
                doubled_grid_shape(3),
                padded_grid_shape(3),
                periodic_z(periodic_z),
                z_period(z_period),
                grid_entire_period(grid_entire_period),
                longitudinal_kicks(longitudinal_kicks),
                commxx_divider_sptr(new Commxx_divider),
                comm2_sptr(),
                comm1_sptr(),
                n_sigma(n_sigma),
                domain_fixed(false),
                have_domains(false),
                diagnostics_list(),
                have_diagnostics(false),
                kick_scale(kick_scale)
{
    constructor_common(grid_shape);
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        Commxx_divider_sptr commxx_divider_sptr,
        std::vector<int > const & grid_shape, bool longitudinal_kicks,
        bool periodic_z, double z_period, bool grid_entire_period,
        double n_sigma, double kick_scale) :
                Collective_operator("space charge 3D open hockney"),
                grid_shape(3),
                doubled_grid_shape(3),
                padded_grid_shape(3),
                periodic_z(periodic_z),
                z_period(z_period),
                grid_entire_period(grid_entire_period),
                longitudinal_kicks(longitudinal_kicks),
                commxx_divider_sptr(commxx_divider_sptr),
                comm2_sptr(),
                comm1_sptr(),
                n_sigma(n_sigma),
                domain_fixed(false),
                have_domains(false),
                diagnostics_list(),
                have_diagnostics(false),
                kick_scale(kick_scale)
{
    constructor_common(grid_shape);
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        Commxx_sptr comm_sptr, std::vector<int > const & grid_shape,
        bool longitudinal_kicks, bool periodic_z, double z_period,
        bool grid_entire_period, double n_sigma, double kick_scale) :
                Collective_operator("space charge 3D open hockney"),
                grid_shape(3),
                doubled_grid_shape(3),
                padded_grid_shape(3),
                periodic_z(periodic_z),
                z_period(z_period),
                grid_entire_period(grid_entire_period),
                longitudinal_kicks(longitudinal_kicks),
                commxx_divider_sptr(new Commxx_divider),
                comm2_sptr(),
                comm1_sptr(),
                n_sigma(n_sigma),
                domain_fixed(false),
                have_domains(false),
                diagnostics_list(),
                have_diagnostics(false),
                kick_scale(kick_scale)
{
    constructor_common(grid_shape);
}

//Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
//        Distributed_fft3d_sptr distributed_fft3d_sptr, bool longitudinal_kicks,
//        bool periodic_z, double z_period, bool grid_entire_period,
//        double n_sigma) :
//        Collective_operator("space charge"), grid_shape(3), doubled_grid_shape(
//                3), padded_grid_shape(3), periodic_z(periodic_z), z_period(
//                z_period), grid_entire_period(grid_entire_period), longitudinal_kicks(
//                longitudinal_kicks), distributed_fft3d_sptr(
//                distributed_fft3d_sptr), comm2_sptr(
//                distributed_fft3d_sptr->get_comm_sptr()), n_sigma(n_sigma), domain_fixed(
//                false), have_domains(false)
//{
//    doubled_grid_shape = distributed_fft3d_sptr->get_shape();
//    for (int i = 0; i < 3; ++i) {
//        grid_shape[i] = doubled_grid_shape[i] / 2;
//    }
//    padded_grid_shape = distributed_fft3d_sptr->get_padded_shape_real();
//    setup_nondoubled_communication();
//    setup_default_options();
//}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney()
{
}

Space_charge_3d_open_hockney *
Space_charge_3d_open_hockney::clone()
{
    return new Space_charge_3d_open_hockney(*this);
}

void
Space_charge_3d_open_hockney::add_diagnostics(Diagnostics_space_charge_3d_hockney_sptr ddiagnostics_sptr)
{
  diagnostics_list.push_back(ddiagnostics_sptr);
  this->have_diagnostics=true;
}

void
Space_charge_3d_open_hockney::set_diagnostics_list(Diagnostics_space_charge_3d_hockneys diagnostics_list)
{
  this->diagnostics_list =diagnostics_list;
  this->have_diagnostics=true;
}

bool
Space_charge_3d_open_hockney::has_diagnostics()
{
    return have_diagnostics;
}

double
Space_charge_3d_open_hockney::get_n_sigma() const
{
    return n_sigma;
}

void
Space_charge_3d_open_hockney::set_green_fn_type(Green_fn_type green_fn_type)
{
    switch (green_fn_type) {
    case pointlike:
        break;
    case linear:
        break;
    default:
        throw runtime_error(
                "Space_charge_3d_open_hockney::set_green_fn_type: invalid green_fn_type");
    }
    this->green_fn_type = green_fn_type;
}

Space_charge_3d_open_hockney::Green_fn_type
Space_charge_3d_open_hockney::get_green_fn_type() const
{
    return green_fn_type;
}

void
Space_charge_3d_open_hockney::set_charge_density_comm(
        Charge_density_comm charge_density_comm)
{
    switch (charge_density_comm) {
    case reduce_scatter:
        break;
    case charge_allreduce:
        break;
    default:
        throw runtime_error(
                "Space_charge_3d_open_hockney::set_charge_density_comm: invalid charge_density_comm");
    }
    this->charge_density_comm = charge_density_comm;
}

Space_charge_3d_open_hockney::Charge_density_comm
Space_charge_3d_open_hockney::get_charge_density_comm() const
{
    return charge_density_comm;
}

void
Space_charge_3d_open_hockney::set_e_field_comm(E_field_comm e_field_comm)
{
    switch (e_field_comm) {
    case gatherv_bcast:
        break;
    case allgatherv:
        break;
    case e_field_allreduce:
        break;
    default:
        throw runtime_error(
                "Space_charge_3d_open_hockney::set_e_field_comm: invalid e_field_comm");
    }

    this->e_field_comm = e_field_comm;
}

Space_charge_3d_open_hockney::E_field_comm
Space_charge_3d_open_hockney::get_e_field_comm() const
{
    return e_field_comm;
}

void
Space_charge_3d_open_hockney::set_doubled_domain()
{
    std::vector<double > doubled_size(3);
    for (int i = 0; i < 3; ++i) {
        doubled_size[i] = 2 * domain_sptr->get_physical_size()[i];
    }
    doubled_domain_sptr = Rectangular_grid_domain_sptr(
            new Rectangular_grid_domain(doubled_size,
                    domain_sptr->get_physical_offset(), doubled_grid_shape,
                    periodic_z));
}

void
Space_charge_3d_open_hockney::set_fixed_domain(
        Rectangular_grid_domain_sptr domain_sptr)
{
    if ((domain_sptr->get_grid_shape()[0] != grid_shape[0])
            || (domain_sptr->get_grid_shape()[1] != grid_shape[1])
            || (domain_sptr->get_grid_shape()[2] != grid_shape[2])) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::set_fixed_domain requires a shape\nequal to that of the parent object, but with zyx ordering.");
    }
    this->domain_sptr = domain_sptr;
    set_doubled_domain();
    domain_fixed = true;
    have_domains = true;
}

// get_smallest_non_tiny is a local function
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

void
Space_charge_3d_open_hockney::update_domain(Bunch const& bunch)
{
    setup_communication(bunch.get_comm_sptr());
    if (!domain_fixed) {
        MArray1d mean(Core_diagnostics::calculate_mean(bunch));
        MArray1d std(Core_diagnostics::calculate_std(bunch, mean));
        std::vector<double > size(3);
        std::vector<double > offset(3);
        const double tiny = 1.0e-10;
        if ((std[Bunch::x] < tiny) && (std[Bunch::y] < tiny)
                && (std[Bunch::z] < tiny)) {
            throw std::runtime_error(
                    "Space_charge_3d_open_hockney::update_domain: all three spatial dimensions have neglible extent");
        }
        if (grid_entire_period) {
            offset[0] = 0.0;
            size[0] = z_period;
        } else {
            offset[0] = mean[Bunch::z];
            size[0] = n_sigma
                    * get_smallest_non_tiny(std[Bunch::z], std[Bunch::x],
                            std[Bunch::y], tiny);
        }
        offset[1] = mean[Bunch::y];
        size[1] = n_sigma
                * get_smallest_non_tiny(std[Bunch::y], std[Bunch::x],
                        std[Bunch::z], tiny);
        offset[2] = mean[Bunch::x];
        size[2] = n_sigma
                * get_smallest_non_tiny(std[Bunch::x], std[Bunch::y],
                        std[Bunch::z], tiny);
        domain_sptr = Rectangular_grid_domain_sptr(
                new Rectangular_grid_domain(size, offset, grid_shape,
                        periodic_z));
        set_doubled_domain();
        have_domains = true;
    }
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_domain_sptr() const
{
    if (!have_domains) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::get_domain_sptr: domain not set");
    }
    return domain_sptr;
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_doubled_domain_sptr() const
{
    if (!have_domains) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::get_doubled_domain_sptr: domain not set");
    }
    return doubled_domain_sptr;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
double t = simple_timer_current();
    update_domain(bunch);
t = simple_timer_show(t, "sc-local-rho-update-domain");
    Rectangular_grid_sptr local_rho_sptr(new Rectangular_grid(domain_sptr));
t = simple_timer_show(t, "sc-local-rho-new");
    deposit_charge_rectangular_zyx(*local_rho_sptr, bunch);
    //deposit_charge_rectangular_zyx_omp_reduce(*local_rho_sptr, bunch);
    //deposit_charge_rectangular_zyx_omp_interleaved(*local_rho_sptr, bunch);
t = simple_timer_show(t, "sc-local-rho-deposit");
    return local_rho_sptr;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_charge_density2_reduce_scatter(
        Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr)
{
    setup_communication(comm_sptr);
// jfa: here is where we do something complicated, but (potentially) efficient
// in calculating a version of the charge density that is just global enough
// to fill in the doubled global charge density

    const double * source = local_charge_density.get_grid_points().origin();

    double * dest;
// dest_array stores the portion of global charge density needed on each
// processor. It has to have the same shape in the non-distributed dimensions
// as the charge density in order to work with MPI_Reduce_scatter.
    MArray3d dest_array(boost::extents[1][1][1]);
    if (real_length > 0) {
        dest_array.resize(
                boost::extents[extent_range(real_lower, real_upper)][grid_shape[1]][grid_shape[2]]);
    }
    dest = multi_array_offset(dest_array, real_lower, 0, 0);
    int error = MPI_Reduce_scatter((void *) source, (void *) dest,
            &real_lengths[0], MPI_DOUBLE, MPI_SUM, comm_sptr->get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney::get_global_charge_density2_reduce_scatter");
    }
    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, doubled_lower,
                    doubled_upper,
                    distributed_fft3d_sptr->get_padded_shape_real(),
                    comm_sptr));
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            for (int k = 0; k < doubled_grid_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] = 0.0;
            }
        }
    }
    for (int i = real_lower; i < real_upper; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] = dest_array[i][j][k];
            }
        }
    }
    return rho2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_charge_density2_allreduce(
        Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr)
{
    setup_communication(comm_sptr);
    int error = MPI_Allreduce(MPI_IN_PLACE,
            (void*) local_charge_density.get_grid_points().origin(),
            local_charge_density.get_grid_points().num_elements(), MPI_DOUBLE,
            MPI_SUM, comm_sptr->get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney::get_global_charge_density2_allreduce");
    }
    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, doubled_lower,
                    doubled_upper,
                    distributed_fft3d_sptr->get_padded_shape_real(),
                    comm_sptr));
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            for (int k = 0; k < doubled_grid_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] = 0.0;
            }
        }
    }
    for (int i = real_lower; i < real_upper; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] =
                        local_charge_density.get_grid_points()[i][j][k];
            }
        }
    }
    return rho2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_charge_density2(
        Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr)
{
    switch (charge_density_comm) {
    case reduce_scatter:
        return get_global_charge_density2_reduce_scatter(local_charge_density,
                comm_sptr);
    case charge_allreduce:
        return get_global_charge_density2_allreduce(local_charge_density,
                comm_sptr);
    default:
        throw runtime_error(
                "Space_charge_3d_open_hockney: invalid charge_density_comm");
    }
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_green_fn2_pointlike()
{
    if (doubled_domain_sptr == NULL) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::get_green_fn2_pointlike called before domain specified");
    }
    int lower = distributed_fft3d_sptr->get_lower();
    int upper = distributed_fft3d_sptr->get_upper();
    Distributed_rectangular_grid_sptr G2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, lower, upper,
                    distributed_fft3d_sptr->get_padded_shape_real(),
                    comm2_sptr));

    double hx = domain_sptr->get_cell_size()[2];
    double hy = domain_sptr->get_cell_size()[1];
    double hz = domain_sptr->get_cell_size()[0];

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
    const double coeff = 2.8;
    double G000 = coeff / std::min(hx, std::min(hy, hz));

    const int num_images = 8;
    int mix, miy; // mirror indices for x- and y-planes
    double dx, dy, dz, G;

// In the following loops we use mirroring for ix and iy, but
// calculate all iz values separately because the mirror points in
// iz may be on another processor.
// Note that the doubling algorithm is not quite symmetric. For
// example, the doubled grid for 4 points in 1d looks like
//    0 1 2 3 4 3 2 1

    #pragma omp parallel for private( dx, dy, dz, G, mix, miy )
    for (int iz = lower; iz < upper; ++iz) {
        if (iz > grid_shape[0]) {
            dz = (doubled_grid_shape[0] - iz) * hz;
        } else {
            dz = iz * hz;
        }
        for (int iy = 0; iy < grid_shape[1] + 1; ++iy) {
            dy = iy * hy;
            miy = doubled_grid_shape[1] - iy;
            if (miy == doubled_grid_shape[1]) {
                miy = iy;
            }
            for (int ix = 0; ix < grid_shape[2] + 1; ++ix) {
                dx = ix * hx;
                mix = doubled_grid_shape[2] - ix;
                if (mix == doubled_grid_shape[2]) {
                    mix = ix;
                }
                if ((ix == 0) && (iy == 0) && (iz == 0)) {
                    G = G000;
                } else {
                    G = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                }
                if (periodic_z) {
                    for (int image = -num_images; image <= num_images;
                            ++image) {
                        if (image != 0) {
                            double dz_image = dz + image * z_period;
                            const double tiny = 1.0e-9;
                            if ((ix == 0) && (iy == 0)
                                    && (std::abs(dz_image) < tiny)) {
                                G += G000;
                            } else {
                                G += 1.0
                                        / sqrt(
                                                dx * dx + dy * dy
                                                        + dz_image * dz_image);
                            }
                        }
                    }
                }
                G2->get_grid_points()[iz][iy][ix] = G;
                // three mirror images
                G2->get_grid_points()[iz][miy][ix] = G;
                G2->get_grid_points()[iz][miy][mix] = G;
                G2->get_grid_points()[iz][iy][mix] = G;
            }
        }
    }

    G2->set_normalization(1.0);

    return G2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_green_fn2_linear()
{
    if (doubled_domain_sptr == NULL) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::get_green_fn2_linear called before domain specified");
    }
    int lower = distributed_fft3d_sptr->get_lower();
    int upper = distributed_fft3d_sptr->get_upper();
    Distributed_rectangular_grid_sptr G2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, lower, upper,
                    distributed_fft3d_sptr->get_padded_shape_real(),
                    comm2_sptr));

    double hx = domain_sptr->get_cell_size()[2];
    double hy = domain_sptr->get_cell_size()[1];
    double hz = domain_sptr->get_cell_size()[0];

    double rr = hx * hx + hy * hy;
    double r1 = sqrt(hx * hx + hy * hy + hz * hz);
    double G000 = (2.0 / rr)
            * (hz * r1 + rr * log((hz + r1) / sqrt(rr)) - hz * hz); // average value of outer cylinder.

    int gz = grid_shape[0];
    int gy = grid_shape[1];
    int gx = grid_shape[2];

    int dgz = doubled_grid_shape[0];
    int dgy = doubled_grid_shape[1];
    int dgx = doubled_grid_shape[2];

    const int num_images = 8;
    int mix, miy; // mirror indices for x- and y-planes
    double x, y, z, G;
    const double epsz = 1.0e-12 * hz;

    #pragma omp parallel for default(none), \
        private( x, y, z, G, mix, miy, rr ), \
        shared( gx, gy, gz, dgx, dgy, dgz, hx, hy, hz, G000, lower, upper, G2 )
    for (int iz = lower; iz < upper; ++iz) {
        z = (iz>gz) ? (dgz-iz)*hz : iz*hz;

        for (int iy = 0; iy <= gy; ++iy) {
            y = iy * hy;
            miy = (iy==gy) ? dgy : (dgy-iy);

            for (int ix = 0; ix <= gx; ++ix) {
                x = ix * hx;
                rr = x * x + y * y;
                mix = (ix==gx) ? dgx : (dgx-ix);

                G = 2.0 * sqrt(rr + z * z) - sqrt(rr + (z - hz) * (z - hz))
                        - sqrt(rr + (z + hz) * (z + hz));
                double T1, T2, r1, r2;
                if (z < -hz) {
                    r1 = (sqrt((z - hz) * (z - hz) + rr) - z + hz)
                            / (sqrt(z * z + rr) - z);
                    T1 = (hz - z) * log(r1);
                    r2 = (sqrt(z * z + rr) - z)
                            / (sqrt((z + hz) * (z + hz) + rr) - z - hz);
                    T2 = (hz + z) * log(r2);
                    G += T1 + T2;
                } else if (std::abs(z + hz) < epsz) {
                    r1 = (sqrt((z - hz) * (z - hz) + rr) - z + hz)
                            / (sqrt(z * z + rr) - z);
                    T1 = (hz - z) * log(r1);
                    G += T1;
                } else if (std::abs(z) < epsz) {
                    if (std::abs(x) + std::abs(y) < 2. * epsz) {
                        G += hz * G000;
                    } /* T1+T2 in fact */else {
                        r1 = (sqrt(hz * hz + rr) + hz) / sqrt(rr);
                        G += 2. * hz * log(r1);
                    }
                } else if (std::abs(z - hz) < epsz) {
                    r1 = (sqrt((z + hz) * (z + hz) + rr) + z + hz)
                            / (sqrt(z * z + rr) + z);
                    T1 = (hz + z) * log(r1);
                    G += T1;
                } else if (z > hz) {
                    r1 = (sqrt(z * z + rr) + z)
                            / (sqrt((z - hz) * (z - hz) + rr) + z - hz);
                    T1 = (hz - z) * log(r1);
                    r2 = (sqrt((z + hz) * (z + hz) + rr) + z + hz)
                            / (sqrt(z * z + rr) + z);
                    T2 = (hz + z) * log(r2);
                    G += T1 + T2;
                } else {
                    throw std::runtime_error(
                            "Space_charge_3d_open_hockney::get_green_fn2 error1");
                }

                if (periodic_z) {
                    throw std::runtime_error(
                            "Space_charge_3d_open_hockney::get_green_fn2_linear: periodic_z not yet implemented");
                    for (int image = -num_images; image < num_images; ++image) {
                        if (image != 0) {
                            double z_image = z + image * z_period;

                            if (z_image < -hz) {
                                r1 = (sqrt((z_image - hz) * (z_image - hz) + rr)
                                        - z_image + hz)
                                        / (sqrt(z_image * z_image + rr)
                                                - z_image);
                                T1 = (hz - z_image) * log(r1);
                                r2 = (sqrt(z_image * z_image + rr) - z_image)
                                        / (sqrt(
                                                (z_image + hz) * (z_image + hz)
                                                        + rr) - z_image - hz);
                                T2 = (hz + z_image) * log(r2);
                                G += T1 + T2;
                            }

                            else if (std::abs(z_image + hz) < epsz) {
                                r1 = (sqrt((z_image - hz) * (z_image - hz) + rr)
                                        - z_image + hz)
                                        / (sqrt(z_image * z_image + rr)
                                                - z_image);
                                T1 = (hz - z_image) * log(r1);
                                G += T1;
                            }

                            else if (std::abs(z_image) < epsz) {
                                if (std::abs(x) + std::abs(y) < 2. * epsz) {
                                    G += hz * G000;
                                } // T1+T2 in fact
                                else {
                                    r1 = (sqrt(hz * hz + rr) + hz) / sqrt(rr);
                                    G += 2. * hz * log(r1);
                                }
                            } else if (std::abs(z_image - hz) < epsz) {
                                r1 = (sqrt((z_image + hz) * (z_image + hz) + rr)
                                        + z_image + hz)
                                        / (sqrt(z_image * z_image + rr)
                                                + z_image);
                                T1 = (hz + z_image) * log(r1);
                                G += T1;
                            } else if (z_image > hz) {
                                r1 = (sqrt(z_image * z_image + rr) + z_image)
                                        / (sqrt(
                                                (z_image - hz) * (z_image - hz)
                                                        + rr) + z_image - hz);
                                T1 = (hz - z_image) * log(r1);
                                r2 = (sqrt((z_image + hz) * (z_image + hz) + rr)
                                        + z_image + hz)
                                        / (sqrt(z_image * z_image + rr)
                                                + z_image);
                                T2 = (hz + z_image) * log(r2);
                                G += T1 + T2;
                            } else {
                                throw std::runtime_error(
                                        "Space_charge_3d_open_hockney::get_green_fn2 error2");
                            }

                        }
                    }
                }

                G2->get_grid_points()[iz][iy][ix] = G;
                // three mirror images
                if (miy < doubled_grid_shape[1]) {
                    G2->get_grid_points()[iz][miy][ix] = G;
                    if (mix < doubled_grid_shape[2]) {
                        G2->get_grid_points()[iz][miy][mix] = G;
                    }
                }
                if (mix < doubled_grid_shape[2]) {
                    G2->get_grid_points()[iz][iy][mix] = G;
                }
            }
        }
    }

    G2->set_normalization(1.0 / (hz * hz));

    return G2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_scalar_field2(
        Distributed_rectangular_grid & charge_density2,
        Distributed_rectangular_grid & green_fn2)
{
    std::vector<int > cshape(
            distributed_fft3d_sptr->get_padded_shape_complex());
    int lower = distributed_fft3d_sptr->get_lower();
    int upper = distributed_fft3d_sptr->get_upper();

    fftw_complex * rho2hat = (fftw_complex*)fftw_malloc(
        sizeof(fftw_complex)*(upper-lower)*cshape[1]*cshape[2]);
    fftw_complex * G2hat   = (fftw_complex*)fftw_malloc(
        sizeof(fftw_complex)*(upper-lower)*cshape[1]*cshape[2]);
    fftw_complex * phi2hat = (fftw_complex*)fftw_malloc(
        sizeof(fftw_complex)*(upper-lower)*cshape[1]*cshape[2]);

    distributed_fft3d_sptr->transform(charge_density2.get_grid_points(), rho2hat);
    distributed_fft3d_sptr->transform(green_fn2.get_grid_points(), G2hat);

    #pragma omp parallel for
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < cshape[1]; ++j) {
            for (int k = 0; k < cshape[2]; ++k) {
                int idx = (i-lower)*cshape[1]*cshape[2] + j*cshape[2] + k;
                phi2hat[idx][0] = rho2hat[idx][0] * G2hat[idx][0] - rho2hat[idx][1] * G2hat[idx][1];
                phi2hat[idx][1] = rho2hat[idx][0] * G2hat[idx][1] + rho2hat[idx][1] * G2hat[idx][0];
            }
        }
    }

    double hx, hy, hz;
    hx = domain_sptr->get_cell_size()[2];
    hy = domain_sptr->get_cell_size()[1];
    hz = domain_sptr->get_cell_size()[0];
    double normalization = hx * hy * hz; // volume element in integral
    normalization *= 1.0 / (4.0 * pi * epsilon0);

    Distributed_rectangular_grid_sptr phi2(
            new Distributed_rectangular_grid(doubled_domain_sptr, lower, upper,
                    distributed_fft3d_sptr->get_padded_shape_real(),
                    comm2_sptr));
    distributed_fft3d_sptr->inv_transform(phi2hat, phi2->get_grid_points());

    normalization *= charge_density2.get_normalization();
    normalization *= green_fn2.get_normalization();
    normalization *= distributed_fft3d_sptr->get_roundtrip_normalization();
    phi2->set_normalization(normalization);

    fftw_free(rho2hat);
    fftw_free(G2hat);
    fftw_free(phi2hat);

    return phi2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::extract_scalar_field(
        Distributed_rectangular_grid const & phi2)
{
    Distributed_rectangular_grid_sptr phi(
            new Distributed_rectangular_grid(domain_sptr, real_doubled_lower,
                    real_doubled_upper, comm1_sptr));

    #pragma omp parallel for
    for (int i = real_doubled_lower; i < real_doubled_upper; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                phi->get_grid_points()[i][j][k] =
                        phi2.get_grid_points()[i][j][k];
            }
        }
    }
    phi->set_normalization(phi2.get_normalization());
    if (comm1_sptr->has_this_rank()) {
        phi->fill_guards();
    }
    return phi;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_electric_field_component(
        Distributed_rectangular_grid const& phi, int component)
{
    int index;
    if (component == 0) {
        index = 2;
    } else if (component == 1) {
        index = 1;
    } else if (component == 2) {
        index = 0;
    } else {
        throw std::runtime_error(
                "Space_charge_3d_open_hockney::get_electric_field_component: component must be 0, 1 or 2");
    }

    Distributed_rectangular_grid_sptr En(
            new Distributed_rectangular_grid(domain_sptr, phi.get_lower(),
                    phi.get_upper(), comm1_sptr));
    MArray3d_ref En_a(En->get_grid_points());
    MArray3d_ref phi_a(phi.get_grid_points());
    int lower_limit, upper_limit;
    if (index == 0) {
        lower_limit = En->get_lower_guard();
        upper_limit = En->get_upper_guard();
    } else {
        lower_limit = 0;
        upper_limit = domain_sptr->get_grid_shape()[index];
    }
    double cell_size = domain_sptr->get_cell_size()[index];
    boost::array<MArray3d::index, 3 > center, left, right;

    #pragma omp parallel for private(center, left, right)
    for (int i = En->get_lower(); i < En->get_upper(); ++i) {
        left[0] = i;
        center[0] = i;
        right[0] = i;
        for (int j = 0; j < domain_sptr->get_grid_shape()[1]; ++j) {
            left[1] = j;
            center[1] = j;
            right[1] = j;
            for (int k = 0; k < domain_sptr->get_grid_shape()[2]; ++k) {
                left[2] = k;
                center[2] = k;
                right[2] = k;

                double delta;
                if (center[index] == lower_limit) {
                    right[index] = center[index] + 1;
                    delta = cell_size;
                } else if (center[index] == upper_limit - 1) {
                    left[index] = center[index] - 1;
                    delta = cell_size;
                } else {
                    right[index] = center[index] + 1;
                    left[index] = center[index] - 1;
                    delta = 2.0 * cell_size;
                }
                // $\vec{E} = - \grad \phi$
                En_a(center) = -(phi_a(right) - phi_a(left)) / delta;
            }
        }
    }
    En->set_normalization(phi.get_normalization());
    return En;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_electric_field_component_gatherv_bcast(
        Distributed_rectangular_grid const& dist_field)
{
    Rectangular_grid_sptr global_field(new Rectangular_grid(domain_sptr));
    const int root = 0;
    int error;
    if (comm1_sptr->has_this_rank()) {
        int rank = comm1_sptr->get_rank();
        error =
                MPI_Gatherv(
                        (void *) (dist_field.get_grid_points().origin()
                                + lowers1[rank]), lengths1[rank], MPI_DOUBLE,
                        (void*) global_field->get_grid_points().origin(),
                        &lengths1[0], &lowers1[0], MPI_DOUBLE, root,
                        comm1_sptr->get());
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Space_charge_3d_open_hockney(MPI_Gatherv)");
        }

    }
    int total_length = grid_shape[0] * grid_shape[1] * grid_shape[2];
    error = MPI_Bcast(global_field->get_grid_points().origin(), total_length,
            MPI_DOUBLE, root, comm2_sptr->get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Bcast)");
    }
    global_field->set_normalization(dist_field.get_normalization());
    return global_field;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_electric_field_component_allgatherv(
        Distributed_rectangular_grid const& dist_field)
{
    Rectangular_grid_sptr global_field(new Rectangular_grid(domain_sptr));
    std::vector<int > lowers12(comm2_sptr->get_size()); // lowers1 on comm2
    std::vector<int > lengths12(comm2_sptr->get_size()); // lengths1 on comm2
    int size1 = lowers1.size();
    for (int rank = 0; rank < comm2_sptr->get_size(); ++rank) {
        if (rank < size1) {
            lowers12[rank] = lowers1[rank];
            lengths12[rank] = lengths1[rank];
        } else {
            lowers12[rank] = 0;
            lengths12[rank] = 0;
        }
    }
    int rank = comm2_sptr->get_rank();
    int error = MPI_Allgatherv(
            (void *) (dist_field.get_grid_points().origin() + lowers12[rank]),
            lengths12[rank], MPI_DOUBLE,
            (void*) global_field->get_grid_points().origin(), &lengths12[0],
            &lowers12[0], MPI_DOUBLE, comm2_sptr->get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Allgatherv)");
    }
    global_field->set_normalization(dist_field.get_normalization());
    return global_field;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_electric_field_component_allreduce(
        Distributed_rectangular_grid const& dist_field)
{
    Rectangular_grid_sptr global_field(new Rectangular_grid(domain_sptr));
    for (int i = 0; i < grid_shape[0]; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                global_field->get_grid_points()[i][j][k] = 0.0;
            }
        }
    }

    std::memset( (void*)global_field->get_grid_points().data(), 0, 
            global_field->get_grid_points().num_elements()*sizeof(double) );

    #pragma omp parallel for
    for (int i = dist_field.get_lower();
            i < std::min(grid_shape[0], dist_field.get_upper()); ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                global_field->get_grid_points()[i][j][k] =
                        dist_field.get_grid_points()[i][j][k];
            }
        }
    }

    int error = MPI_Allreduce(MPI_IN_PLACE,
            (void*) global_field->get_grid_points().origin(),
            global_field->get_grid_points().num_elements(), MPI_DOUBLE, MPI_SUM,
            comm2_sptr->get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Allreduce in get_global_electric_field_component_allreduce)");
    }
    global_field->set_normalization(dist_field.get_normalization());
    return global_field;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_electric_field_component(
        Distributed_rectangular_grid const& dist_field)
{
    Rectangular_grid_sptr global_electric_field_sptr;
    if (e_field_comm==gatherv_bcast) {
        global_electric_field_sptr=get_global_electric_field_component_gatherv_bcast(dist_field);
    }
    else if(e_field_comm==allgatherv){
        global_electric_field_sptr=get_global_electric_field_component_allgatherv(dist_field);
    }
    else if(e_field_comm==e_field_allreduce){
        global_electric_field_sptr=get_global_electric_field_component_allreduce(dist_field);
    }
    else throw runtime_error(
                 "Space_charge_3d_open_hockney: invalid e_field_comm");


//AM!  make sure the field is zero at the edge of the grid
// THIS toghether with zero charge distribution at the edge of the grid is essential for a conservative approximation
    MArray3d_ref grid_points(global_electric_field_sptr->get_grid_points());
    if (!global_electric_field_sptr->get_domain_sptr()->is_periodic()){
          for (int k=0; k<grid_points.shape()[2];++k){
              for (int j=0; j<grid_points.shape()[1];++j){
                  grid_points[0][j][k]=0.;
                  grid_points[grid_points.shape()[0]-1][j][k]=0.;      
              }
          } 
    }
    for (int i=0; i<grid_points.shape()[0];++i){
        for (int j=0; j<grid_points.shape()[1];++j){
            grid_points[i][j][0]=0.;
            grid_points[i][j][grid_points.shape()[2]-1]=0.;
            
        }
    } 
   
    for (int i=0; i<grid_points.shape()[0];++i){
        for (int k=0; k<grid_points.shape()[2];++k){
            grid_points[i][0][k]=0.;         
            grid_points[i][grid_points.shape()[1]-1][k]=0.;
        }
    }
    return global_electric_field_sptr; 
}

void
Space_charge_3d_open_hockney::set_kick_scale(double ks)
{
    kick_scale = ks;
}

double
Space_charge_3d_open_hockney::get_kick_scale() const
{
    return kick_scale;
}

void
Space_charge_3d_open_hockney::do_diagnostics(Rectangular_grid const& En, int component, double time_step, Step & step, 
                                          Bunch & bunch)
{   
   if (have_diagnostics) {
      if ((component==0) || (component==1)){
         double step_beta=step.get_betas()[component];
         for (Diagnostics_space_charge_3d_hockneys::const_iterator d_it = diagnostics_list.begin();
            d_it != diagnostics_list.end(); ++d_it){   
            if (bunch.is_bucket_index_assigned()){
                if ((*d_it)->get_bunch().get_bucket_index()==bunch.get_bucket_index()){                  
                    (*d_it)->update(bunch, En, component, time_step, step_beta); 
                    if (component==1) (*d_it)->write();
                }
            }
            else{
                    (*d_it)->update(bunch, En, component, time_step, step_beta); 
                    if (component==1) (*d_it)->write();
            }                                                 
         }
      }    
   } 
   
}  


void
Space_charge_3d_open_hockney::apply_kick(Bunch & bunch,
        Rectangular_grid const& En, double delta_t, int component)
{
  
  
  //AM: kicks  in the z_lab frame 
 //Delta p_x&=& F_x \Delta t&=& - q \frac{1}{\gamma^2} \frac{\partial \Phi'}{\partial x} \Delta t=q \frac{1}{\beta \gamma^2} E_{grid~x} \Delta t\\
 //Delta E &= & q E_z \Delta s&=& q \frac{1}{\gamma^2 \beta} \frac{\partial \Phi'}{\partial ct} \beta c\Delta t=-q \frac{c}{\beta \gamma^2 }E_{grid~z} \Delta t\\
 // 1/beta factor in E_grid from charge deposition on (x,y,cdt) coordinates grid,  \Phi'=\frac{1}{\beta}\Phi_{grid}

 
 

    bunch.convert_to_state(Bunch::fixed_z_lab);
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double gamma=bunch.get_reference_particle().get_gamma();
    double beta=bunch.get_reference_particle().get_beta();
// unit_conversion: [kg m/s] to [Gev/c]
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
// scaled p = p/p_ref
    double p_ref=bunch.get_reference_particle().get_momentum();
    double factor = kick_scale * unit_conversion * q * delta_t* En.get_normalization()/
            (p_ref*gamma*gamma*beta); // transverse kicks
   

    int ps_component = 2 * component + 1;
    Rectangular_grid_domain & domain(*En.get_domain_sptr());
    MArray3d_ref grid_points(En.get_grid_points());

    if (component==2)
    {
        factor *= -p_ref; 
        double m = bunch.get_mass();

        #pragma omp parallel for
        for (int part = 0; part < bunch.get_local_num(); ++part) 
        {
            double x = bunch.get_local_particles()[part][Bunch::x];
            double y = bunch.get_local_particles()[part][Bunch::y];
            double z = bunch.get_local_particles()[part][Bunch::z];   

            double grid_val = interpolate_rectangular_zyx(x, y, z, domain, grid_points);

            double p = p_ref +bunch.get_local_particles()[part][Bunch::dpop] * p_ref;        
            double Eoc_i = std::sqrt(p * p + m * m);
            double Eoc_f = Eoc_i + factor * grid_val;
            double delta_dpop = (std::sqrt(Eoc_f*Eoc_f-m*m) - std::sqrt(Eoc_i*Eoc_i-m*m))/p_ref;

            bunch.get_local_particles()[part][ps_component] += delta_dpop;
        }      

        // spectator particles
        #pragma omp parallel for
        for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
        {
            double x = bunch.get_local_spectator_particles()[part][Bunch::x];
            double y = bunch.get_local_spectator_particles()[part][Bunch::y];
            double z = bunch.get_local_spectator_particles()[part][Bunch::z];   

            double grid_val = interpolate_rectangular_zyx(x, y, z, domain, grid_points);

            double p = p_ref + bunch.get_local_spectator_particles()[part][Bunch::dpop] * p_ref;        
            double Eoc_i = std::sqrt(p * p + m * m);
            double Eoc_f = Eoc_i + factor * grid_val;
            double delta_dpop = (std::sqrt(Eoc_f*Eoc_f-m*m) - std::sqrt(Eoc_i*Eoc_i-m*m))/p_ref;

            bunch.get_local_spectator_particles()[part][ps_component] += delta_dpop;
        }      
    }
    else
    {
        #pragma omp parallel for
        for (int part = 0; part < bunch.get_local_num(); ++part) 
        {
              double x = bunch.get_local_particles()[part][Bunch::x];
              double y = bunch.get_local_particles()[part][Bunch::y];
              double z = bunch.get_local_particles()[part][Bunch::z];   

              double grid_val = interpolate_rectangular_zyx(x, y, z, domain, grid_points);      

              bunch.get_local_particles()[part][ps_component] += factor * grid_val;
        }

        // spectator particles
        #pragma omp parallel for
        for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
        {
              double x = bunch.get_local_spectator_particles()[part][Bunch::x];
              double y = bunch.get_local_spectator_particles()[part][Bunch::y];
              double z = bunch.get_local_spectator_particles()[part][Bunch::z];   

              double grid_val = interpolate_rectangular_zyx(x, y, z, domain, grid_points);      

              bunch.get_local_spectator_particles()[part][ps_component] += factor * grid_val;
        }
    }
}

void
Space_charge_3d_open_hockney::apply(Bunch & bunch, double time_step,
        Step & step, int verbosity, Logger & logger)
{
    if (bunch.get_total_num() > 1) {
        double t = simple_timer_current();
        setup_communication(bunch.get_comm_sptr());
        int comm_compare;
        MPI_Comm_compare(comm2_sptr->get(), bunch.get_comm().get(),
                &comm_compare);
        if ((comm_compare == MPI_UNEQUAL)
                && (charge_density_comm != charge_allreduce)) {
            throw std::runtime_error(
                    "Space_charge_3d_open_hockney: set_charge_density_comm(charge_allreduce) required when comm != bunch comm");
        }
        t = simple_timer_show(t, "sc-setup-communication");

     
        bunch.convert_to_state(Bunch::fixed_z_lab);
        
        t = simple_timer_show(t, "sc-convert-to-state");
        Rectangular_grid_sptr local_rho(get_local_charge_density(bunch)); // [C/m^3]
        t = simple_timer_show(t, "sc-get-local-rho");
        Distributed_rectangular_grid_sptr rho2(
                get_global_charge_density2(*local_rho, bunch.get_comm_sptr())); // [C/m^3]
        t = simple_timer_show(t, "sc-get-global-rho");
        local_rho.reset();
        Distributed_rectangular_grid_sptr G2; // [1/m]
        if (green_fn_type == pointlike) {
            G2 = get_green_fn2_pointlike();
        } else if (green_fn_type == linear) {
            G2 = get_green_fn2_linear();
        } else {
            throw std::runtime_error(
                    "Space_charge_3d_open_hockney::apply: unknown green_fn_type");
        }
        t = simple_timer_show(t, "sc-get-green-fn");
        Distributed_rectangular_grid_sptr phi2(get_scalar_field2(*rho2, *G2)); // [V]
        t = simple_timer_show(t, "sc-get-phi2");
        rho2.reset();
        G2.reset();
        Distributed_rectangular_grid_sptr phi(extract_scalar_field(*phi2));
        t = simple_timer_show(t, "sc-get-phi");
      //  bunch.periodic_sort(Bunch::z); // A.M this causes troubles
        t = simple_timer_show(t, "sc-sort");
        phi2.reset();
        int max_component;
        if (longitudinal_kicks) {
            max_component = 3;
        } else {
            max_component = 2;
        }
        for (int component = 0; component < max_component; ++component) {
            Distributed_rectangular_grid_sptr local_En(
                    get_electric_field_component(*phi, component)); // [V/m]
            t = simple_timer_show(t, "sc-get-local-en");
            Rectangular_grid_sptr En(
                    get_global_electric_field_component(*local_En)); // [V/m]
            t = simple_timer_show(t, "sc-get-global-en");
            do_diagnostics(*En,component, time_step,step, bunch);
            apply_kick(bunch, *En, time_step, component);
            t = simple_timer_show(t, "sc-apply-kick");
        }
    }
}
#endif

