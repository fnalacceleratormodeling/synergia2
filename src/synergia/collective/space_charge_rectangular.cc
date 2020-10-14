#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/core_diagnostics.h"
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
                        //<< hgrid(z*gx*gy + y*gx + x) << ", ";
                        << hgrid(x*gy*gz + y*gz + z) << ", ";
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

    struct alg_phi
    {
        karray1d_dev phi;
        int lower, gy, gz;
        double igygz, igz;
        double px, py, pz;
        double gamma;

        alg_phi(karray1d_dev const& phi,
                int lower, int gy, int gz,
                double px, double py, double pz,
                double gamma)
            : phi(phi)
            , lower(lower), gy(gy), gz(gz)
            , igygz(1.0/(gy*gz))
            , igz(1.0/gz)
            , px(px), py(py), pz(pz)
            , gamma(gamma)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            const int ix = i * igygz;
            const int iy = (i - ix*gy*gz) * igz;
            const int iz = i - ix*gy*gz - iy*gz;

            int xt = ix + lower + 1;
            int yt = iy + 1;

            double denominator 
                = mconstants::pi * mconstants::pi * ( 
                        xt*xt / (px*px) + 
                        yt*yt / (py*py) + 
                        4.0 * iz*iz / (pz*pz*gamma*gamma) );

            int base = ix*gy*gz + iy*gz + iz;

            phi(base*2+0) /= denominator;
            phi(base*2+1) /= denominator;

        }
    };

    struct alg_force_extractor
    {
        karray1d_dev phi;
        karray1d_dev enx;
        karray1d_dev eny;
        karray1d_dev enz;

        int gx, gy, gz;
        int dgx, dgy;
        double ihx, ihy, ihz;
        double igygz;
        double igz;

        alg_force_extractor(
                karray1d_dev const& phi,
                karray1d_dev const& enx,
                karray1d_dev const& eny,
                karray1d_dev const& enz,
                std::array<int, 3> const& g,
                std::array<double, 3> const& h )
            : phi(phi), enx(enx), eny(eny), enz(enz)
            , gx(g[0]), gy(g[1]), gz(g[2])
            , ihx(0.5/h[0]), ihy(0.5/h[1]), ihz(0.5/h[2])
            , igygz(1.0/(gy*gz))
            , igz(1.0/gz)
        { }


        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            int ix = i * igygz;
            int iy = (i - ix*gy*gz) * igz;
            int iz = i - ix*gy*gz - iy*gz;

            int ixl, ixr, iyl, iyr, izl, izr;
            double idx, idy, idz;

            // all x-y boundaries will be skipped (set to 0)
            if (ix==0 || ix==gx-1 || iy==0 || iy==gy-1)
                return;

            ixl = ix - 1; ixr = ix + 1;
            iyl = iy - 1; iyr = iy + 1;
            izl = iz - 1; izr = iz + 1;

            idx = ihx; idy = ihy; idz = ihz;

            // periodic z
            if (iz==0) { izl = gz-1; }
            else if (iz==gz-1) { izr = 0; }

            int idx_r = ixr*gy*gz + iy*gz + iz;
            int idx_l = ixl*gy*gz + iy*gz + iz;
            enx(i) = -(phi(idx_r) - phi(idx_l)) * idx;

            idx_r = ix*gy*gz + iyr*gz + iz;
            idx_l = ix*gy*gz + iyl*gz + iz;
            eny(i) = -(phi(idx_r) - phi(idx_l)) * idy;

            idx_r = ix*gy*gz + iy*gz + izr;
            idx_l = ix*gy*gz + iy*gz + izl;
            enz(i) = -(phi(idx_r) - phi(idx_l)) * idz;
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

                if (ix>=0 && ix<gx-1 && iy>=0 && iy<gy-1)
                {
                    while (iz > gz-1) iz -= gz;
                    while (iz < 0) iz += gz;

                    int izp1 = (iz==(gz-1)) ? 0 : iz + 1;

                    double val = 0;
                    int base = 0;

                    // enz
                    base = ix*gy*gz + iy*gz;
                    val  = aox * aoy * aoz * enz(base+iz);      // x, y, z
                    val += aox * aoy *  oz * enz(base+izp1);    // x, y, z+1
                    val += aox *  oy * aoz * enz(base+gz+iz);   // x, y+1, z
                    val += aox *  oy *  oz * enz(base+gz+izp1); // x, y+1, z+1

                    base = (ix+1)*gy*gz + iy*gz;
                    val += ox * aoy * aoz * enz(base+iz);       // x+1, y, z
                    val += ox * aoy *  oz * enz(base+izp1);     // x+1, y, z+1
                    val += ox *  oy * aoz * enz(base+gz+iz);    // x+1, y+1, z
                    val += ox *  oy *  oz * enz(base+gz+izp1);  // x+1, y+1, z+1

                    double p = pref + parts(i, 5) * pref;
                    double Eoc_i = std::sqrt(p*p+m*m);
                    double Eoc_f = Eoc_i - factor * pref * val;
                    double dpop = (std::sqrt(Eoc_f*Eoc_f-m*m) - std::sqrt(Eoc_i*Eoc_i-m*m))/pref;

                    parts(i, 5) += dpop;

                    // eny
                    base = ix*gy*gz + iy*gz;
                    val  = aox * aoy * aoz * eny(base+iz);      // x, y, z
                    val += aox * aoy *  oz * eny(base+izp1);    // x, y, z+1
                    val += aox *  oy * aoz * eny(base+gz+iz);   // x, y+1, z
                    val += aox *  oy *  oz * eny(base+gz+izp1); // x, y+1, z+1

                    base = (ix+1)*gy*gz + iy*gz;
                    val += ox * aoy * aoz * eny(base+iz);       // x+1, y, z
                    val += ox * aoy *  oz * eny(base+izp1);     // x+1, y, z+1
                    val += ox *  oy * aoz * eny(base+gz+iz);    // x+1, y+1, z
                    val += ox *  oy *  oz * eny(base+gz+izp1);  // x+1, y+1, z+1

                    parts(i, 3) += factor * val;

                    // enx
                    base = ix*gy*gz + iy*gz;
                    val  = aox * aoy * aoz * enx(base+iz);      // x, y, z
                    val += aox * aoy *  oz * enx(base+izp1);    // x, y, z+1
                    val += aox *  oy * aoz * enx(base+gz+iz);   // x, y+1, z
                    val += aox *  oy *  oz * enx(base+gz+izp1); // x, y+1, z+1

                    base = (ix+1)*gy*gz + iy*gz;
                    val += ox * aoy * aoz * enx(base+iz);       // x+1, y, z
                    val += ox * aoy *  oz * enx(base+izp1);     // x+1, y, z+1
                    val += ox *  oy * aoz * enx(base+gz+iz);    // x+1, y+1, z
                    val += ox *  oy *  oz * enx(base+gz+izp1);  // x+1, y+1, z+1

                    parts(i, 1) += factor * val;
                }
            }
        }
    };

}




Space_charge_rectangular::Space_charge_rectangular(
            Space_charge_rectangular_options const& ops)
    : Collective_operator("sc_rectangular", 1.0)
    , options(ops)
    , domain(ops.shape, {1.0, 1.0, 1.0})
    , ffts()
{
}


void 
Space_charge_rectangular::apply_impl(
            Bunch_simulator& sim, 
            double time_step, 
            Logger& logger)
{
    logger << "    Space charge 3d open hockney\n";

    scoped_simple_timer timer("sc_rect_total");

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
Space_charge_rectangular::apply_bunch(
            Bunch& bunch, 
            Distributed_fft3d_rect& fft,
            double time_step, 
            Logger& logger)
{
    update_domain(bunch);

    get_local_charge_density(bunch);
    get_global_charge_density(bunch);

    double gamma = bunch.get_reference_particle().get_gamma();

    get_local_phi(fft, gamma);
    get_global_phi(fft);

    auto fn_norm = get_normalization_force();

    extract_force();
    apply_kick(bunch, fn_norm, time_step);
}

void
Space_charge_rectangular::construct_workspaces(
        Bunch_simulator const& sim)
{
    // shape
    auto const& s = options.shape;

    // fft objects
    for(size_t t=0; t<2; ++t)
    {
        int num_local_bunches = sim[t].get_bunch_array_size();
        ffts[t] = std::vector<Distributed_fft3d_rect>(num_local_bunches);

        for(size_t b=0; b<num_local_bunches; ++b)
        {
            auto comm = sim[t][b]
                .get_comm()
                .divide(options.comm_group_size);

            ffts[t][b].construct(s, comm);
        }
    }

    // local workspaces
    int nz_cplx = Distributed_fft3d_rect::get_padded_shape_cplx(s[2]);

    rho = karray1d_dev("rho", s[0] * s[1] * s[2]);
    phi = karray1d_dev("phi", s[0] * s[1] * s[2]);

    phihat = karray1d_dev("phihat", s[0] * s[1] * nz_cplx*2);

    h_rho = Kokkos::create_mirror_view(rho);
    h_phi = Kokkos::create_mirror_view(phi);

    // En is in the original domain
    enx = karray1d_dev("enx", s[0] * s[1] * s[2]);
    eny = karray1d_dev("eny", s[0] * s[1] * s[2]);
    enz = karray1d_dev("enz", s[0] * s[1] * s[2]);
}


void
Space_charge_rectangular::update_domain(Bunch const& bunch)
{
    double beta = bunch.get_reference_particle().get_beta();
    auto dsize = options.pipe_size;
    dsize[2] /= beta;    // size in z_lab frame, longitudinal cdt coordinate

    // A.M physical_offsets of the domain should be rescaled too, 
    // but in this case they are zero    
    domain = Rectangular_grid_domain(
            options.shape, dsize, {0.0, 0.0, 0.0}, true);
}

void
Space_charge_rectangular::get_local_charge_density(Bunch const& bunch)
{
    scoped_simple_timer timer("sc_rect_local_rho");

    auto g = domain.get_grid_shape();
    //g[2] = Distributed_fft3d::get_padded_shape_real(g[2]);
    //g[2] = (g[2]/2+1)*2;

#ifdef Kokkos_ENABLE_OPENMP
    deposit_charge_rectangular_3d_omp_reduce_xyz(rho,
            domain, g, bunch);
#else
    deposit_charge_rectangular_3d_kokkos_scatter_view_xyz(rho,
            domain, g, bunch);
#endif
}


void
Space_charge_rectangular::get_global_charge_density(Bunch const& bunch)
{
    // do nothing if the bunch occupis a single rank
    if (bunch.get_comm().size() == 1) return;

    scoped_simple_timer timer("sc_rect_global_rho");

    auto g = domain.get_grid_shape();

    simple_timer_start("sc_rect_global_rho_copy");
    Kokkos::deep_copy(h_rho, rho);
    simple_timer_stop("sc_rect_global_rho_copy");

    simple_timer_start("sc_rect_global_rho_reduce");
    int err = MPI_Allreduce( MPI_IN_PLACE,
                             (void*)h_rho.data(), 
                             h_rho.extent(0),
                             MPI_DOUBLE, 
                             MPI_SUM, 
                             bunch.get_comm() );
    simple_timer_stop("sc_rect_global_rho_reduce");

    if (err != MPI_SUCCESS)
    {
        throw std::runtime_error( 
                "MPI error in Space_charge_rectangular"
                "(MPI_Allreduce in get_global_charge_density)" );
    }

    simple_timer_start("sc_rect_global_rho_copy");
    Kokkos::deep_copy(rho, h_rho);
    simple_timer_stop("sc_rect_global_rho_copy");
}


void
Space_charge_rectangular::get_local_phi(
        Distributed_fft3d_rect& fft, double gamma)
{
    scoped_simple_timer timer("sc_rect_local_phi");

    alg_zeroer az{phihat};
    Kokkos::parallel_for(phihat.extent(0), az);

    fft.transform(rho, phihat);

    int lower = fft.get_lower();
    int upper = fft.get_upper();
    auto ps = options.pipe_size;
    auto g = domain.get_grid_shape();
    int gy = g[1];
    int gz_padded_cplx = g[2]/2 + 1;

    alg_phi aphi(phihat, lower, gy, gz_padded_cplx,
            ps[0], ps[1], ps[2], gamma);

    Kokkos::parallel_for(
            (upper-lower)*gy*gz_padded_cplx, aphi);

    fft.inv_transform(phihat, phi);
}

void
Space_charge_rectangular::get_global_phi(
        Distributed_fft3d_rect const& fft)
{
    // do nothing if the solver only has a single rank
    if (fft.get_comm().size() == 1) return;

    scoped_simple_timer timer("sc_rect_global_phi");

    Kokkos::deep_copy(h_phi, phi);

    int err = MPI_Allreduce( MPI_IN_PLACE,
                             (void*)h_phi.data(), 
                             h_phi.extent(0),
                             MPI_DOUBLE, 
                             MPI_SUM, 
                             fft.get_comm() );

    if (err != MPI_SUCCESS)
    {
        throw std::runtime_error( 
                "MPI error in Space_charge_3d_open_hockney"
                "(MPI_Allreduce in get_global_electric_force2_allreduce)" );
    }

    Kokkos::deep_copy(phi, h_phi);
}

double
Space_charge_rectangular::get_normalization_force()
{
    auto g = domain.get_grid_shape();
    return 1.0 / (4.0 * g[0] * g[1] * g[2] * pconstants::epsilon0);
}

void
Space_charge_rectangular::extract_force()
{
    scoped_simple_timer timer("sc_rect_get_en");

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();

    // phi is in (gx, gy, gz)
    // en{x|y|z} is in (gx, gy, gz)
    alg_force_extractor alg(phi, enx, eny, enz, g, h);
    Kokkos::parallel_for(g[0]*g[1]*g[2], alg);
    Kokkos::fence();
}

void
Space_charge_rectangular::apply_kick(
        Bunch& bunch,
        double fn_norm,
        double time_step)
{
    scoped_simple_timer timer("sc_rect_kick");

    auto ref = bunch.get_reference_particle();

    double q = bunch.get_particle_charge() * pconstants::e;
    double m = bunch.get_mass();

    double gamma = ref.get_gamma();
    double beta  = ref.get_beta();
    double pref  = ref.get_momentum();

    double unit_conversion = pconstants::c / (1e9 * pconstants::e);
    double factor = unit_conversion * q * time_step * fn_norm 
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

