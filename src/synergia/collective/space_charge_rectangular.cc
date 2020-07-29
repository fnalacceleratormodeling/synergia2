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
    , fft()
    , comm(Commxx::Null)
{
}


void 
Space_charge_rectangular::apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger)
{
    logger << "    Space charge 3d open hockney\n";

    scoped_simple_timer timer("sc3d_total");
    apply_bunch(simulator[0][0], time_step, logger);
}

void
Space_charge_rectangular::apply_bunch(
            Bunch& bunch, 
            double time_step, 
            Logger& logger)
{
    set_domain(bunch);

    //Logger l;
    //print_statistics(bunch, l);

    get_local_charge_density(bunch);
    get_global_charge_density(bunch);

    double gamma = bunch.get_reference_particle().get_gamma();

    get_local_phi(gamma);
    get_global_phi();

    //Logger l2;
    //print_grid(l2, phi, 32, 36, 32, 33, 32, 33, 64, 64, 130);

    auto fn_norm = get_normalization_force();

    extract_force();
    apply_kick(bunch, fn_norm, time_step);

    //Logger l3;
    //print_statistics(bunch, l3);
}

void
Space_charge_rectangular::construct_workspaces(
        std::array<int, 3> const& s)
{
    int nz_cplx = fft.padded_nz_cplx();

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
Space_charge_rectangular::set_domain(Bunch const& bunch)
{
    if (comm.is_null())
    {
        comm = bunch.get_comm().divide(options.comm_group_size);
        fft.construct(options.shape, comm);
        construct_workspaces(options.shape);
    }

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
    deposit_charge_rectangular_3d_kokkos_scatter_view(rho,
            domain, g, bunch);
#endif

#if 0
    Rectangular_grid_sptr rho_sptr(new Rectangular_grid(domain_sptr));
    deposit_charge_rectangular_xyz(*rho_sptr, bunch);

    int error = MPI_Allreduce(MPI_IN_PLACE, 
            (void*)  rho_sptr->get_grid_points().origin(), 
            rho_sptr->get_grid_points().num_elements(), 
            MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
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
Space_charge_rectangular::get_local_phi(double gamma)
{
    alg_zeroer az{phihat};
    Kokkos::parallel_for(phihat.extent(0), az);

    //Logger l;
    //print_grid(l, rho, 32, 36, 32, 33, 32, 33, 64, 64, 130);

    fft.transform(rho, phihat);
    //print_grid(l, phihat, 32, 36, 32, 33, 32, 33, 64, 64, 130);

    int lower = fft.get_lower();
    int upper = fft.get_upper();

    auto ps = options.pipe_size;

    auto g = domain.get_grid_shape();
    int gy = g[1];
    int gz_padded_cplx = g[2]/2 + 1;

    for(int x=lower; x<upper; ++x)
    {
        int xt = x + 1;

        for(int y=0; y<gy; ++y)
        {
            int yt = y + 1;

            for(int z=0; z<gz_padded_cplx; ++z)
            {
                double denominator = pi * pi * (
                    xt * xt / (ps[0] * ps[0]) +
                    yt * yt / (ps[1] * ps[1]) +
                    4.0 * z * z / (ps[2] * ps[2] * gamma * gamma)
                );

                int base = x*gy*gz_padded_cplx + y*gz_padded_cplx + z;

                phihat(base*2+0) /= denominator;
                phihat(base*2+1) /= denominator;
            }
        }
    }

    fft.inv_transform(phihat, phi);
    //phi_local->set_normalization(1./(4.*shape[0]*shape[1]*shape[2]*epsilon0));
}

void
Space_charge_rectangular::get_global_phi()
{
    // do nothing if the solver only has a single rank
    if (comm.size() == 1) return;

    scoped_simple_timer timer("sc_rect_global_phi");

    Kokkos::deep_copy(h_phi, phi);

    int err = MPI_Allreduce( MPI_IN_PLACE,
                             (void*)h_phi.data(), 
                             h_phi.extent(0),
                             MPI_DOUBLE, 
                             MPI_SUM, 
                             comm );

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

#if 0
    Logger l;
    print_grid(l, enx, 32, 36, 32, 33, 32, 33, 64, 64, 130);
    print_grid(l, eny, 32, 36, 32, 33, 32, 33, 64, 64, 130);
    print_grid(l, enz, 32, 36, 32, 33, 32, 33, 64, 64, 130);
#endif
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




#if 0
Space_charge_rectangular::Space_charge_rectangular(Commxx_sptr comm_f_sptr, std::vector<double > const & pipe_size, 
			std::vector<int > const & grid_shape, bool equally_spread):
Collective_operator("space_charge_rectangular"), pipe_size(pipe_size), 
grid_shape(grid_shape),  comm_f_sptr(comm_f_sptr), equally_spread(equally_spread), have_domain(false),
have_diagnostics(false)//,fftw_helper_sptr(),domain_sptr()
{

 try{  
    this->have_fftw_helper=false;
    construct_fftw_helper(comm_f_sptr);
     if ((!comm_f_sptr->has_this_rank()) && (equally_spread)) throw std::runtime_error(
		  "Space_charge_rectangular:: equally_spread is incompatible with this choice of comm_f_sptr ");

 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}


Space_charge_rectangular::Space_charge_rectangular(std::vector<double > const & pipe_size, std::vector<int > const & grid_shape):
Collective_operator("space_charge_rectangular"),  pipe_size(pipe_size),  grid_shape(grid_shape), 
have_domain(false),have_diagnostics(false)//,comm_f_sptr(),fftw_helper_sptr(),domain_sptr()
{

 try{   
    this->have_fftw_helper=false;
    this->equally_spread=false; 

 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}
 
Space_charge_rectangular::Space_charge_rectangular()
{
}

Space_charge_rectangular *
Space_charge_rectangular::clone()
{
    return new Space_charge_rectangular(*this);
}


template<class Archive>
    void
    Space_charge_rectangular::save(Archive & ar, const unsigned int version) const
    {
       
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(comm_f_sptr)
	    &  BOOST_SERIALIZATION_NVP(grid_shape)
	    &  BOOST_SERIALIZATION_NVP(pipe_size) 
	    &  BOOST_SERIALIZATION_NVP(have_fftw_helper)
	    &  BOOST_SERIALIZATION_NVP(equally_spread)
        &  BOOST_SERIALIZATION_NVP(have_diagnostics)
        &  BOOST_SERIALIZATION_NVP(diagnostics_list);
    }

template<class Archive>
    void
    Space_charge_rectangular::load(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
        ar & BOOST_SERIALIZATION_NVP(comm_f_sptr)
        &  BOOST_SERIALIZATION_NVP(grid_shape)
        &  BOOST_SERIALIZATION_NVP(pipe_size) 
        &  BOOST_SERIALIZATION_NVP(have_fftw_helper)
        &  BOOST_SERIALIZATION_NVP(equally_spread)
        &  BOOST_SERIALIZATION_NVP(have_diagnostics)
        &  BOOST_SERIALIZATION_NVP(diagnostics_list);
        domain_sptr.reset();
        have_domain=false;       
        if (have_fftw_helper) { this->have_fftw_helper=false;          
                                construct_fftw_helper(comm_f_sptr);  
                               }
                                       
    }

template
void
Space_charge_rectangular::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Space_charge_rectangular::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Space_charge_rectangular::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Space_charge_rectangular::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);





BOOST_CLASS_EXPORT_IMPLEMENT(Space_charge_rectangular)


Space_charge_rectangular::~Space_charge_rectangular()
{
 //fftw_mpi_cleanup();
}

void
Space_charge_rectangular::construct_fftw_helper(Commxx_sptr comm_sptr)
{
      if (!have_fftw_helper){
	 if (comm_sptr->has_this_rank()){
	    this->fftw_helper_sptr=  Fftw_rectangular_helper_sptr (  new   Fftw_rectangular_helper(grid_shape, comm_sptr));
	    this->have_fftw_helper=true;
	 }
	  this->comm_f_sptr=comm_sptr;
      }
      else {
	  throw std::runtime_error(
		  "Space_charge_rectangular::construct_fftw_helper:   already has fftw_helper ");
      }
    
}


void
Space_charge_rectangular::set_fftw_helper(Commxx_sptr comm_sptr, bool equally_spread)
{

 try{
    if (!have_fftw_helper){
	construct_fftw_helper(comm_sptr);
    }
    else {
	if (comm_sptr->has_this_rank()) this->fftw_helper_sptr->reset_comm_f(comm_sptr);
	this->comm_f_sptr=comm_sptr;	    
   }
   this->equally_spread=equally_spread; 
   if ((!comm_sptr->has_this_rank()) && (equally_spread)) throw std::runtime_error(
		  "Space_charge_rectangular:: set fftw: equally_spread is incompatible with this choice of comm_sptr ");
 }
 catch (std::exception const& e){
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 111);
    }
}

bool  
Space_charge_rectangular::get_have_fftw_helper() const
{
  return have_fftw_helper;
}


Commxx_sptr 
Space_charge_rectangular::get_comm_sptr() const
{
  return comm_f_sptr;
}  

std::vector<double >
Space_charge_rectangular::get_pipe_size() const
{
return pipe_size;
}

std::vector<int >
Space_charge_rectangular::get_grid_shape() const
{
return grid_shape;
}

void
Space_charge_rectangular::add_diagnostics(Diagnostics_space_charge_rectangular_sptr ddiagnostics_sptr)
{
  diagnostics_list.push_back(ddiagnostics_sptr);
  this->have_diagnostics=true;
}

void
Space_charge_rectangular::set_diagnostics_list(Diagnostics_space_charge_rectangulars diagnostics_list)
{
  this->diagnostics_list =diagnostics_list;
  this->have_diagnostics=true;
}



Rectangular_grid_domain_sptr
Space_charge_rectangular::get_domain_sptr() const
{
return domain_sptr;
}

void
Space_charge_rectangular::set_domain(Bunch const & bunch)
{
  if (!have_domain){
        double beta = bunch.get_reference_particle().get_beta();
        std::vector<double >  dsize=pipe_size;
        dsize[2] /= beta;    // size in z_lab frame, longitudinal cdt coordinate
         // A.M physical_offsets of the domain should be rescaled too, but in this case they are zero    
        this->domain_sptr = Rectangular_grid_domain_sptr(
                      new Rectangular_grid_domain(dsize,  grid_shape, true));
         have_domain=true;  
  }
}

bool 
Space_charge_rectangular::get_have_domain() const
{
  return have_domain;
}

Fftw_rectangular_helper_sptr
Space_charge_rectangular::get_fftw_helper_sptr() const
{
    return fftw_helper_sptr;
}



Rectangular_grid_sptr
Space_charge_rectangular::get_charge_density(Bunch const& bunch)
{
    double t;
    t = simple_timer_current(); 

    if (!have_domain) set_domain(bunch);   
    Rectangular_grid_sptr rho_sptr(new Rectangular_grid(domain_sptr));
    deposit_charge_rectangular_xyz(*rho_sptr, bunch);
    //t = simple_timer_show(t, "sc_get_charge_density: depozit_xyz");

    t = simple_timer_current();
    int error = MPI_Allreduce(MPI_IN_PLACE, (void*)  rho_sptr->get_grid_points().origin(),
                               rho_sptr->get_grid_points().num_elements(), MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());


     t = simple_timer_show(t, "sc_get_charge_density: allmpireduce");

    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_rectangular: MPI_Allreduce in get_charge_density");
    }


    return rho_sptr;
}



Distributed_rectangular_grid_sptr
Space_charge_rectangular::get_phi_local( Rectangular_grid & rho, double const& gamma)
{

    if (!have_fftw_helper)  throw std::runtime_error(
                "Space_charge_rectangular::get_phi_local  space_charge does not have have_fftw_helper defined");

    if (!comm_f_sptr->has_this_rank()) throw 
             std::runtime_error("space charge rectangular, get_phi_local, comm_f_sptr has no rank");
//    double t;
//    t = simple_timer_current();

    MArray3d_ref rho_ref(rho.get_grid_points());
    std::vector<int > shape(domain_sptr->get_grid_shape());

    ptrdiff_t local_nx=get_fftw_helper_sptr()->get_local_nx();
    ptrdiff_t local_x_start=get_fftw_helper_sptr()->get_local_x_start();



    int lower=local_x_start;
    int upper=local_x_start+local_nx;
    std::string solver("rectangular");
    std::vector<double > local_physical_size(domain_sptr->get_physical_size());
    std::vector<double > local_physical_offset(domain_sptr->get_physical_offset());
    // local_physical_size[0]=domain_sptr->get_cell_size()[0]*local_nx;
    //local_physical_offset[0]=domain_sptr->get_left()[0]+local_x_start*domain_sptr->get_cell_size()[0];
    std::vector<int > shape_local(shape);
    shape_local[0]=local_nx;
    Distributed_rectangular_grid_sptr phi_local(
	new Distributed_rectangular_grid(local_physical_size, local_physical_offset, shape_local,
	    true, lower, upper, comm_f_sptr, solver)
	    );


    MArray3d_ref phi_local_ref(phi_local->get_grid_points());
    MArray3d_ref rho_local_ref(rho_ref.origin()+local_x_start*shape[1]*shape[2],boost::extents[local_nx][shape[1]][shape[2]]);

    // t = simple_timer_current();
    get_fftw_helper_sptr()->transform(rho_local_ref, phi_local_ref);

  // t = simple_timer_show(t, "sc_get_phi_local: fftw_dst_direct");


    

    const int memory_fudge_factor = 1;
    fftw_complex *rho_nmp_local;
    rho_nmp_local= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *local_nx*shape[1]*(shape[2]/2+1)*memory_fudge_factor);
    int dim[] = {shape[2]};
    fftw_plan plan=fftw_plan_many_dft_r2c(1, dim, local_nx*shape[1], phi_local_ref.origin(), NULL,
				      1, shape[2], rho_nmp_local, NULL, 1,shape[2]/2+1,
					FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
  // t = simple_timer_show(t, "sc_get_phi_local: fftw1d_direct");


    MArray3dc_ref rho_nmp_ref(reinterpret_cast<std::complex<double>*>(rho_nmp_local), boost::extents[local_nx][shape[1]][shape[2]/2+1]);

    for (int n=0; n < local_nx; ++n){
	int nt=n+1+local_x_start;
	for (int m=0; m< shape[1]; ++m){
	    int mt=m+1;
	    for (int p=0; p< shape[2]/2+1; ++p){
		double denominator=pi*pi*
		    (nt*nt/(pipe_size[0]*pipe_size[0])+mt*mt/(pipe_size[1]*pipe_size[1])+4.*p*p/(pipe_size[2]*pipe_size[2]*gamma*gamma));
			rho_nmp_ref[n][m][p] /= denominator; // delta Phi =- rho
	    }
	}
    }
    //t = simple_timer_show(t, "sc_get_phi_local: loop_phi_nmp");

    plan=fftw_plan_many_dft_c2r(1, dim, local_nx*shape[1], rho_nmp_local, NULL, 1,shape[2]/2+1,
					phi_local_ref.origin(), NULL, 1,shape[2],
					FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(rho_nmp_local);


    //t = simple_timer_show(t, "sc_get_phi_local: fftw1d_inverse");
    //t = simple_timer_current();
    get_fftw_helper_sptr()->inv_transform(phi_local_ref,phi_local_ref);
  // t = simple_timer_show(t, "sc_get_phi_local: fftw_dst_inverse");

    phi_local->set_normalization(1./(4.*shape[0]*shape[1]*shape[2]*epsilon0));
    return phi_local;
 }



void
Space_charge_rectangular::fill_guards_pplanes(Distributed_rectangular_grid & phi, int lower, int upper, int lengthx,
                          MArray2d & g_lower, MArray2d &g_upper)
{
  
   if (!comm_f_sptr->has_this_rank()) throw 
             std::runtime_error("space charge rectangular, fill_guards_pplanes, comm_f_sptr has no rank");

    int mpi_compare;
    MPI_Comm_compare(phi.get_comm().get(), comm_f_sptr->get(), &mpi_compare) ;
    if  (mpi_compare != MPI_IDENT)    {
        throw std::runtime_error("space charge rectangular, fill_guards_pplanes, phi comm is not the same as space_charge comm_f");
    }

    int lrank=comm_f_sptr->get_rank();
    std::vector<int > shape_phi(phi.get_domain().get_grid_shape());
    int message_size = shape_phi[1] * shape_phi[2];
    int shapex=upper-lower;

    void *send_buffer, *recv_buffer;
    MPI_Status status;

     // send to the right

    if ((upper < lengthx) &&  (upper >0)) {
        send_buffer=reinterpret_cast<void*>(phi.get_grid_points().origin()+(shapex-1)*message_size);
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, lrank + 1, lrank, comm_f_sptr->get());
    }
    if (lower > 0) {
        recv_buffer=reinterpret_cast<void*>(g_lower.data());
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, lrank - 1, lrank - 1,
                 comm_f_sptr->get(), &status);
    }
 // send to the left

    if (lower > 0) {
        send_buffer=reinterpret_cast<void*>(phi.get_grid_points().origin());
        MPI_Send(send_buffer, message_size, MPI_DOUBLE, lrank - 1, lrank,
                 comm_f_sptr->get());
    }
    if ((upper < lengthx) &&  (upper >0)){
        recv_buffer=reinterpret_cast<void*>(g_upper.data());
        MPI_Recv(recv_buffer, message_size, MPI_DOUBLE, lrank + 1, lrank + 1,
                 comm_f_sptr->get(), &status);
    }
}


Rectangular_grid_sptr
Space_charge_rectangular::get_En(Distributed_rectangular_grid &phi, int component)
{
  
    if (!comm_f_sptr->has_this_rank()) throw 
             std::runtime_error("space charge rectangular, get_En, comm_f_sptr has no rank");

if ((component < 0) || (component > 2)) {
        std::stringstream message("");
        message << "calculate_E_n: invalid argument component=" << component
                << ". Argument be in range 0<=component<=2";
        throw std::invalid_argument(message.str());
    }
    double t;
    t = simple_timer_current();


    std::vector<int > shape_phi(phi.get_domain().get_grid_shape());
    MArray3d E_local(boost::extents[shape_phi[0]][shape_phi[1]][shape_phi[2]]);
    std::vector<double > hi(domain_sptr->get_cell_size());
    double h(hi[component]);


    int size=comm_f_sptr->get_size();
    int lrank=comm_f_sptr->get_rank();
    std::vector<int > shape(domain_sptr->get_grid_shape());

    int lower =  phi.get_lower();
    int upper =  phi.get_upper();
    int lengthx=shape[0];

    MArray2d guard_lower(boost::extents[shape[1]][shape[2]]);
    MArray2d guard_upper(boost::extents[shape[1]][shape[2]]);


    Rectangular_grid_domain_sptr domain_local_sptr(phi.get_domain_sptr());
    Rectangular_grid_sptr En_local(new Rectangular_grid(domain_local_sptr));
    MArray3d_ref En_local_a(En_local->get_grid_points());
    MArray3d_ref phi_a(phi.get_grid_points());
    boost::array<MArray3d::index, 3 > center, left, right;


    if (component==0) {
        fill_guards_pplanes(phi,  lower, upper, lengthx, guard_lower, guard_upper);
        for (int i = 0; i < upper-lower; ++i) {
                left[0] = i;
                center[0] = i;
                right[0] = i;
                for (int j = 0; j < shape[1]; ++j) {
                    left[1] = j;
                    center[1] =j;
                    right[1] = j;
                    for (int k = 0; k < shape[2]; ++k) {
                        left[2] = k;
                        center[2] = k;
                        right[2] = k;
                        right[component] = std::min(int(center[component] + 1), upper-lower-1);
                        left[component] = std::max(int(center[component] - 1),0);
                        double delta=  (right[component]-left[component])*h;
                        double phi_right, phi_left;
                        if ((center[0]==upper-lower-1) && (upper < shape[0])) {
                            phi_right=guard_upper[j][k];
                            delta=2*h;
                        }
                        else  phi_right=phi_a(right);
                        if ((center[0]==0) && (lower>0)){
                           phi_left=guard_lower[j][k];
                           delta=2*h;
                        }
                        else phi_left=phi_a(left);
                         //  $\vec{E} = - \grad \phi$
                        En_local_a(center)=-(phi_right - phi_left) / delta;
                     }
                 }
            }


    }
    else if (component==2) { // we have periodicity
        for (int i = 0; i < upper-lower; ++i) {
            left[0] = i;
            center[0] = i;
            right[0] = i;
            for (int j = 0; j < shape[1]; ++j) {
                left[1] = j;
                center[1] =j;
                right[1] = j;
                for (int k = 0; k < shape[2]; ++k) {
                    left[2] = k;
                    center[2] = k;
                    right[2] = k;
                    right[component] = (center[component] + 1) % shape[2];
                    left[component] = (center[component]-1 + shape[2])% shape[2];
                    double delta=2.*h;
                      //  $\vec{E} = - \grad \phi$
                    En_local_a(center)= -(phi_a(right) - phi_a(left)) / delta;
                 }
             }
        }
    }
    else {
        for (int i = 0; i < upper-lower; ++i) {
            left[0] = i;
            center[0] = i;
            right[0] = i;
            for (int j = 0; j < shape[1]; ++j) {
                left[1] = j;
                center[1] =j;
                right[1] = j;
                for (int k = 0; k < shape[2]; ++k) {
                    left[2] = k;
                    center[2] = k;
                    right[2] = k;
                    right[component] = std::min(int(center[component] + 1),shape[component] - 1);
                    left[component] = std::max(int(center[component] - 1),0);
                    double delta=(right[component]-left[component])*h;
                      //  $\vec{E} = - \grad \phi$
                    En_local_a(center)= -(phi_a(right) - phi_a(left)) / delta;
                 }
             }
        }
    }

  //  t = simple_timer_show(t, "get_En: calculate E local");

    Rectangular_grid_sptr En(new Rectangular_grid(domain_sptr));

     std::vector<int> uppers(phi.get_uppers());
     std::vector<int> receive_counts(phi.get_lengths()), receive_offsets(size);
     for (int i=0; i< size; ++i) {
         receive_offsets.at(i) = uppers.at(i)*shape[1]*shape[2]-receive_counts.at(i);
     }

     t = simple_timer_current();

    if (equally_spread){ 
	int error = MPI_Allgatherv(reinterpret_cast<void*>(En_local_a.origin()),
		  receive_counts[lrank], MPI_DOUBLE,
		  reinterpret_cast<void*>(En->get_grid_points().origin()),
					  &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, comm_f_sptr->get());

	if (error != MPI_SUCCESS) {
	    throw std::runtime_error(
		"MPI error in Space_charge_rectangular(MPI_Allgatherv in get_En: En_local)");
	}
	En->set_normalization(phi.get_normalization()); // we should have here  \div $\vec{E}=rho/epsilon   
    }
    else{
	int error = MPI_Gatherv(reinterpret_cast<void*>(En_local_a.origin()),
		  receive_counts[lrank], MPI_DOUBLE,
		  reinterpret_cast<void*>(En->get_grid_points().origin()),
					  &receive_counts[0], &receive_offsets[0], MPI_DOUBLE, 0,comm_f_sptr->get());

	if (error != MPI_SUCCESS) {
	    throw std::runtime_error(
		"MPI error in Space_charge_rectangular(MPI_Gatherv in get_En: En_local)");
	}
	
    }   
    
   //AM!  make sure the field is zero at the edge of the grid
// THIS toghether with zero charge distribution at the edge of the grid is essential for a conservative approximation
    MArray3d_ref grid_points(En->get_grid_points());
    for (int j=0; j<grid_points.shape()[1];++j){
        for (int k=0; k<grid_points.shape()[2];++k){
            grid_points[0][j][k]=0.;
            grid_points[grid_points.shape()[0]-1][j][k]=0.;
            
        }
    }    
    for (int i=0; i<grid_points.shape()[0];++i){
        for (int k=0; k<grid_points.shape()[2];++k){
            grid_points[i][0][k]=0.;         
            grid_points[i][grid_points.shape()[1]-1][k]=0.;
        }
    } 
    
    
    
    t = simple_timer_show(t, "get_En:  gather En");
     
    
    
    
    return En;

}

void
Space_charge_rectangular::do_diagnostics(Rectangular_grid const& En, int component, double time_step, Step & step, 
                                          Bunch & bunch)
{   
   if (have_diagnostics) {
      if ((component==0) || (component==1)){
         double step_beta=step.get_betas()[component];
         for (Diagnostics_space_charge_rectangulars::const_iterator d_it = diagnostics_list.begin();
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
Space_charge_rectangular::apply_kick(Bunch & bunch, Rectangular_grid const& En, double  delta_t, int component)
{
  
 //AM: kicks  in the z_lab frame 
 //Delta p_x&=& F_x \Delta t&=& - q \frac{1}{\gamma^2} \frac{\partial \Phi'}{\partial x} \Delta t=q \frac{1}{\beta \gamma^2} E_{grid~x} \Delta t\\
 //Delta E &= & q E_z \Delta s&=& q \frac{1}{\gamma^2 \beta} \frac{\partial \Phi'}{\partial ct} \beta c\Delta t=-q \frac{c}{\beta \gamma^2 }E_{grid~z} \Delta t\\
 // 1/beta factor in E_grid from charge deposition on (x,y,cdt) coordinates grid 

 
 

    bunch.convert_to_state(Bunch::fixed_z_lab);
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double gamma=bunch.get_reference_particle().get_gamma();
    double beta=bunch.get_reference_particle().get_beta();
// unit_conversion: [kg m/s] to [Gev/c] 
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
// scaled p = p/p_ref
    double p_ref=bunch.get_reference_particle().get_momentum();
    double factor = unit_conversion * q * delta_t* En.get_normalization()/
            (p_ref*gamma*gamma*beta); // transverse kicks
   
      
    int ps_component = 2 * component + 1;
  
 
    Rectangular_grid_domain & domain(*En.get_domain_sptr());
    MArray3d_ref grid_points(En.get_grid_points());
   
    if (component==2){
       factor *= -p_ref; 
       double m = bunch.get_mass();
       for (int part = 0; part < bunch.get_local_num(); ++part) {
          double x = bunch.get_local_particles()[part][Bunch::x];
          double y = bunch.get_local_particles()[part][Bunch::y];
          double z = bunch.get_local_particles()[part][Bunch::z];   
          double grid_val =  interpolate_rectangular_xyz(x, y, z, domain,
                grid_points);
          double p=p_ref +bunch.get_local_particles()[part][Bunch::dpop] * p_ref;        
          double Eoc_i = std::sqrt(p * p + m * m);
          double Eoc_f=  Eoc_i + factor * grid_val;
          double delta_dpop=(std::sqrt(Eoc_f*Eoc_f-m*m)-std::sqrt(Eoc_i*Eoc_i-m*m))/p_ref;
          bunch.get_local_particles()[part][ps_component] += delta_dpop;
       }      
    }
    else{
        for (int part = 0; part < bunch.get_local_num(); ++part) {
          double x = bunch.get_local_particles()[part][Bunch::x];
          double y = bunch.get_local_particles()[part][Bunch::y];
          double z = bunch.get_local_particles()[part][Bunch::z];        
          double grid_val = interpolate_rectangular_xyz(x, y, z, domain,
                  grid_points);   
          bunch.get_local_particles()[part][ps_component] += factor * grid_val; 
        }
    }
   
}

std::vector<Rectangular_grid_sptr>
Space_charge_rectangular::get_Efield(Rectangular_grid & rho,Bunch const& bunch, int max_component, double const & gamma )
{	
   if (equally_spread) throw std::runtime_error
              	("Space_charge_rectangular get_Efield: don't call this function for true equally_spread ");
  
   std::vector<Rectangular_grid_sptr> Efield;
    if (comm_f_sptr->has_this_rank()){
       Distributed_rectangular_grid_sptr phi_local(get_phi_local(rho,gamma)); // \nabla phi= -rho/epsilon0; [phi]=kg*m^2*C^{-1}*s^{-2}            
       for (int component = 0; component < max_component; ++component) {
          Efield.push_back(get_En(*phi_local, component));
       }	
    } 
    else{
      for (int component = 0; component < max_component; ++component){ 
        Efield.push_back(Rectangular_grid_sptr(new Rectangular_grid(domain_sptr)));			
      }
      int mpi_compare;
      MPI_Comm_compare(comm_f_sptr->get_parent_sptr()->get(), bunch.get_comm_sptr()->get(), &mpi_compare);
      if ((mpi_compare != MPI_IDENT) && ( mpi_compare != MPI_CONGRUENT)){	       
                  throw std::runtime_error
                  ("Space_charge_rectangular get_Efield: comm_f_sptr parent and bunch.comm are not congruent");
      } 
    }  // comm_f_sptr->has_this_rank()
   
   
   
   

   // cast Efield from rank=0 of comm_spc to whole bunch.get_comm
   std::vector<int > shape(domain_sptr->get_grid_shape());
   int count=shape[0]*shape[1]*shape[2];
   for (int component = 0; component < max_component; ++component){
      int error=MPI_Bcast(Efield[component]->get_grid_points().origin(), count, 
				  MPI_DOUBLE, 0, bunch.get_comm().get());
	if (error != MPI_SUCCESS) {
		throw std::runtime_error(
		  "MPI error in Space_charge_rectangular, get_Efield: MPI_Bcast Efield failed)");
	}	  
	Efield[component]->set_normalization(1./(4.*shape[0]*shape[1]*shape[2]*epsilon0));
   }    
   return Efield;   
}  


void
Space_charge_rectangular::apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger)
{ 
    double t,t1;
   // t = simple_timer_current();   
    bunch.convert_to_state(Bunch::fixed_z_lab);
    set_domain(bunch);
    Rectangular_grid_sptr rho_sptr(get_charge_density(bunch)); // [C/m^3], [Q/DxDyD(ct)] in z_lab frame
    t = simple_timer_show(t, "sc_apply: get-rho");
    double gamma=bunch.get_reference_particle().get_gamma();
    int max_component(3);   
    if (equally_spread){     
        if (comm_f_sptr->get_parent_sptr().get()!= 0) {
	    int mpi_compare;
	    MPI_Comm_compare(comm_f_sptr->get_parent_sptr()->get(), bunch.get_comm_sptr()->get(), &mpi_compare);
            if ((mpi_compare != MPI_IDENT) && ( mpi_compare != MPI_CONGRUENT)){
               throw std::runtime_error
               ("Space_charge_rectangular apply, equally_spread=1: comm_f_sptr parent and bunch.comm are not congruent");
            } 	  
	    }	   
        Distributed_rectangular_grid_sptr phi_local(get_phi_local(*rho_sptr,gamma)); // Phi_local is Phi'(\gamma z)=\Phi(z), but the grid is (x,y,cdt)
	    for (int component = 0; component < max_component; ++component) {	  
	        Rectangular_grid_sptr  En(get_En(*phi_local, component)); // E=-/grad phi; [E]=kg*m/(C*s^2)=N/C	
            do_diagnostics(*En,component, time_step,step, bunch);
	        apply_kick(bunch, *En, time_step, component);	 
	   }
    }
    else{
        std::vector<Rectangular_grid_sptr>  Efield(get_Efield(*rho_sptr, bunch, max_component,gamma));
        for (int component = 0; component < max_component; ++component) { 
            do_diagnostics(*(Efield[component]),component, time_step,step, bunch);
            apply_kick(bunch, *(Efield[component]), time_step, component);  
        }
    }      
    
     t = simple_timer_show(t, "sc_apply: 3x apply_kick and get En");
     t1 = simple_timer_show(t1, "sc_aplly total");
}
#endif

