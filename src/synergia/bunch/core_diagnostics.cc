#include "core_diagnostics.h"
#include <cmath>
//#include "Eigen/Core"
//#include "Eigen/LU"
#include <stdexcept>
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/synergia_omp.h"

#include <functional>

//using namespace Eigen;

#if 0
// note: cannot get the template F to work with CUDA
//
template<typename F, size_t... I>
struct particle_reducer
{
    typedef double value_type[];
    typedef Particles::size_type size_type;

    constexpr size_type value_count;
    F f;
    ConstParticles p;

    particle_reducer(F f, ConstParticles const & parts)
        : value_count(sizeof...(I)), f(f), p(parts)
    { }

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type i, value_type sum) const
    {
        //sum[0] += f(p, i, 0);
        //sum[2] += f(p, i, 2);
        //sum[4] += f(p, i, 4);

        // it works, pros? cons?
        //(void)(std::initializer_list<double>{ (sum[I] += f(p, i, I))... } );

        // only with c++17
        //((sum[I] += f(p, i, I)),...);

        constexpr std::array<size_t, sizeof...(I)> indices{{I...}};
        for (size_type j=0; j<value_count; ++j) sum[j] += f(p, i, indices[j]);
    }

    KOKKOS_INLINE_FUNCTION void
    join(volatile value_type dst, const volatile value_type src) const
    {
        for (size_type j=0; j<value_count; ++j) dst[j] += src[j];
    }

    KOKKOS_INLINE_FUNCTION void
    init(value_type sum) const
    {
        for (size_type j=0; j<value_count; ++j) sum[j] = 0.0;
    }
};

using p_fun_t = double(ConstParticles const &, int, int);
using fun_t = std::function<double(ConstParticles const &, int, int)>;
#endif

namespace core_diagnostics_impl
{
    struct mean_tag   { constexpr static int size = 6; };
    struct z_mean_tag { constexpr static int size = 1; };
    struct std_tag    { constexpr static int size = 6; };
    struct mom2_tag   { constexpr static int size = 36; };

    template<typename F>
    struct particle_reducer
    {
        typedef double value_type[];

        const int value_count = F::size;
        ConstParticles p;
        karray1d_dev dev_mean;

        particle_reducer(
                ConstParticles const & p, 
                karray1d const & mean = karray1d("mean", 6) )
            : p(p), dev_mean("dev_mean", 6) 
        { Kokkos::deep_copy(dev_mean, mean); }

#if 0
        KOKKOS_INLINE_FUNCTION
        void join(volatile value_type dst, const volatile value_type src) const
        { for (int j=0; j<value_count; ++j) dst[j] += src[j]; }
#endif

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, value_type sum) const;
    };

    template<>
    KOKKOS_INLINE_FUNCTION
    void particle_reducer<mean_tag>::operator()   (const int i, value_type sum) const
    { 
        for (int j=0; j<value_count; ++j) sum[j] += p(i, j); 
    }

    template<>
    KOKKOS_INLINE_FUNCTION
    void particle_reducer<z_mean_tag>::operator() (const int i, value_type sum) const
    { 
        sum[0] += p(i, 4);
    }

    template<>
    KOKKOS_INLINE_FUNCTION
    void particle_reducer<std_tag>::operator()    (const int i, value_type sum) const
    { 
        for (int j=0; j<value_count; ++j) 
            sum[j] += (p(i, j) - dev_mean(j)) * (p(i, j) - dev_mean(j)); 
    }

    template<>
    KOKKOS_INLINE_FUNCTION
    void particle_reducer<mom2_tag>::operator()   (const int i, value_type sum) const
    {
        for (int j=0; j<6; ++j)
        {
            double diff_j = p(i, j) - dev_mean(j);
            for (int k=0; k<=j; ++k)
            {
                double diff_k = p(i, k) - dev_mean(k);
                sum[j*6+k] += diff_j * diff_k;
            }
        }
    }
}

karray1d
Core_diagnostics::calculate_mean(Bunch const & bunch)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::mean_tag;

    karray1d mean("mean", 6);

    auto particles = bunch.get_local_particles();
    const int npart = bunch.get_local_num();

    particle_reducer<mean_tag> pr(particles);
    Kokkos::parallel_reduce("cal_mean", npart, pr, mean.data());

    Kokkos::fence();

    double t = simple_timer_current();
    MPI_Allreduce(MPI_IN_PLACE, mean.data(), 6, MPI_DOUBLE, MPI_SUM, bunch.get_comm());
    t = simple_timer_show(t, "allmpireduce_in_diagnostic mean");

    for (int i=0; i<6; ++i) mean(i) /= bunch.get_total_num();

    return mean;
}

double
Core_diagnostics::calculate_z_mean(Bunch const& bunch)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::z_mean_tag;

    double mean = 0;

    auto particles = bunch.get_local_particles();
    const int npart = bunch.get_local_num();

    particle_reducer<z_mean_tag> pr(particles);
    Kokkos::parallel_reduce(npart, pr, &mean);

    Kokkos::fence();

    double t = simple_timer_current();
    MPI_Allreduce(MPI_IN_PLACE, &mean, 1, MPI_DOUBLE, MPI_SUM, bunch.get_comm());
    t = simple_timer_show(t, "allmpireduce_in_diagnostic z mean");

    mean = mean / bunch.get_total_num();

    return mean;
}

double
Core_diagnostics::calculate_z_std(Bunch const& bunch, double const& mean)
{
    double sum = 0;
    double std = 0;
#if 0
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double diff = particles[part][4] - mean;
        sum += diff * diff;
    }
    MPI_Allreduce(&sum, &std, 1, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    std = std::sqrt(std / bunch.get_total_num());
#endif
    return std;
}

karray1d
Core_diagnostics::calculate_spatial_mean(Bunch const& bunch)
{
    karray1d mean("mean", 3);
#if 0
    double sum[3] = { 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    int npart = bunch.get_local_num();

    #pragma omp parallel shared(npart, particles)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int l = npart / nt;
        int s = it * l;
        int e = (it==nt-1) ? npart : (it+1)*l;

        double lsum0 = 0, lsum1 = 0, lsum2 = 0;

        for (int part = s; part < e; ++part)
        {
            lsum0 += particles[part][0];
            lsum1 += particles[part][2];
            lsum2 += particles[part][4];
        }

        #pragma omp critical
        {
            sum[0] += lsum0;
            sum[1] += lsum1;
            sum[2] += lsum2;
        }
    }

    double t;
    t = simple_timer_current();
    MPI_Allreduce(sum, mean.origin(), 3, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    t = simple_timer_show(t, "allmpireduce_in_diagnostic mean");

    for (int i = 0; i < 3; ++i) {
        mean[i] /= bunch.get_total_num();
    }
#endif
    return mean;
}

karray1d
Core_diagnostics::calculate_std(Bunch const & bunch, karray1d const & mean)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::std_tag;

    karray1d std("std", 6);

    auto particles = bunch.get_local_particles();
    const int npart = bunch.get_local_num();

    particle_reducer<std_tag> pr(particles, mean);
    Kokkos::parallel_reduce(npart, pr, std.data());

    double t = simple_timer_current();
    MPI_Allreduce(MPI_IN_PLACE, std.data(), 6, MPI_DOUBLE, MPI_SUM, bunch.get_comm());
    t = simple_timer_show(t, "allmpireduce_in_diagnostic mean");

    for (int i=0; i<6; ++i) std(i) = std::sqrt(std(i) / bunch.get_total_num());

    return std;
}

karray1d
Core_diagnostics::calculate_spatial_std(Bunch const& bunch, karray1d const& mean)
{
    karray1d std("std", 3);
#if 0
    double sum[3] = { 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    int npart = bunch.get_local_num();

    #pragma omp parallel shared(npart, particles)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int l = npart / nt;
        int s = it * l;
        int e = (it==nt-1) ? npart : (it+1)*l;

        double lsum[3] = { 0, 0, 0 };
        double diff;

        for(int part = s; part < e; ++part)
        {
            diff = particles[part][0] - mean[0]; lsum[0] += diff * diff;
            diff = particles[part][2] - mean[1]; lsum[1] += diff * diff;
            diff = particles[part][4] - mean[2]; lsum[2] += diff * diff;
        }
 
        #pragma omp critical
        {
            sum[0] += lsum[0];
            sum[1] += lsum[1];
            sum[2] += lsum[2];
        } // end of omp critical
    }

    MPI_Allreduce(sum, std.origin(), 3, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 3; ++i) {
        std[i] = std::sqrt(std[i] / bunch.get_total_num());
    }
#endif
    return std;
}

karray2d
Core_diagnostics::calculate_mom2(Bunch const & bunch, karray1d const & mean)
{
    using core_diagnostics_impl::particle_reducer;
    using core_diagnostics_impl::mom2_tag;

    karray2d mom2("mom2", 6, 6);
    karray1d sum2("sum2", 36);

    auto particles = bunch.get_local_particles();
    auto npart = bunch.get_local_num();

    particle_reducer<mom2_tag> pr(particles, mean);
    Kokkos::parallel_reduce(npart, pr, sum2.data());

    for (int i=0; i<5; ++i)
        for (int j=i+1; j<6; ++j)
            sum2(i*6+j) = sum2(j*6+i);

    MPI_Allreduce(MPI_IN_PLACE, sum2.data(), 36, MPI_DOUBLE, MPI_SUM, bunch.get_comm());

    for (int i=0; i<6; ++i)
        for (int j=0; j<6; ++j)
            mom2(i, j) = sum2(i*6+j) / bunch.get_total_num();

    return mom2;

#if 0
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            sum2[i][j] = 0.0;
        }
    }
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            double diff_i = particles[part][i] - mean[i];
            for (int j = 0; j <= i; ++j) {
                double diff_j = particles[part][j] - mean[j];
                sum2[i][j] += diff_i * diff_j;
            }
        }
    }
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            sum2[i][j] = sum2[j][i];
        }
    }
    MPI_Allreduce(sum2.origin(), mom2.origin(), 36, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            mom2[i][j] = mom2[j][i] = mom2[i][j] / bunch.get_total_num();
        }
    }
#endif

}

karray1d
Core_diagnostics::calculate_min(Bunch const& bunch)
{
    karray1d min("min", 3);
#if 0
    double lmin[3] = { 1.0e100, 1.0e100, 1.0e100 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        if (particles[part][0] < lmin[0]) {
            lmin[0] = particles[part][0];
        }
        if (particles[part][2] < lmin[1]) {
            lmin[1] = particles[part][2];
        }
        if (particles[part][4] < lmin[2]) {
            lmin[2] = particles[part][4];
        }

    }
    MPI_Allreduce(lmin, min.origin(), 3, MPI_DOUBLE, MPI_MIN,
            bunch.get_comm().get());

#endif
    return min;
}

karray1d
Core_diagnostics::calculate_max(Bunch const& bunch)
{
    karray1d max("max", 3);
#if 0
    double lmax[3] = { -1.0e100, -1.0e100, -1.0e100 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        if (particles[part][0] > lmax[0]) {
            lmax[0] = particles[part][0];
        }
        if (particles[part][2] > lmax[1]) {
            lmax[1] = particles[part][2];
        }
        if (particles[part][4] > lmax[2]) {
            lmax[2] = particles[part][4];
        }

    }
    MPI_Allreduce(lmax, max.origin(), 3, MPI_DOUBLE, MPI_MAX,
            bunch.get_comm().get());

#endif
    return max;
}

void
Core_diagnostics::print_bunch_parameters(karray2d const& mom2, double beta)
{
#if 0
  ///emitx,emity,emitz correspond to sigma^2/beta for a matched beam. Note there is no pi factor in our definition. 
  /// 95% emitts...corespond to Fermilab measured emittances defined as (6 pi sigma^2/beta0.

   double gamma=1./sqrt(1.-beta*beta);
   double pz = gamma * beta *pconstants::mp;
   double energy=pconstants::mp * gamma;
   std::vector<double> units(6);
   units[0]=1.;
   units[1]=1./pz;
   units[2]=1.;
   units[3]=1./pz;
   units[4]=1./beta;
   units[5]=1./pz;
   
   double emitx=sqrt(mom2[0][0]*mom2[1][1]-mom2[0][1]*mom2[1][0])/units[0]/units[1]; // this is xrms^2/beta_lattice !!!!!
   double emity=sqrt(mom2[2][2]*mom2[3][3]-mom2[2][3]*mom2[3][2])/units[2]/units[3]; // this is yrms^2/beta_lattice !!!!!
   double emitz=sqrt(mom2[4][4]*mom2[5][5]-mom2[4][5]*mom2[5][4])/units[4]/units[5]; // this is zrms^2/beta_lattice !!!!!

   // don't print out 10000 copies of the the bunch parameters
   Logger logger(0);

  logger<< "************ BEAM MATCHED PARAMETERS*************"<<std::endl;
  logger<<"*    emitx="<< emitx<<" meters*GeV/c   ="<<emitx/pz<<" meters*rad (synergia units)="<<emitx/pz/mconstants::pi<<" pi*meters*rad"<<std::endl;
  logger<<"*    emity="<< emity<< " meters*GeV/c   ="<<emity/pz<< " meters*rad (synergia units)="<< emity/pz/mconstants::pi<<" pi*meters*rad"<<std::endl;
  logger<<"*    emitz="<<emitz<<" meters*GeV/c ="<<emitz*1.e9/(pconstants::c)<<" eV*s =" <<emitz*beta*beta*energy/pz
            <<" meters*GeV ="<<emitz/pz/beta<<" [cdt*dp/p] (synergia units)"<<std::endl;
  logger<<std::endl;
  logger<<"*    90%emitx="<< 4.605*mconstants::pi*emitx/pz<<"  meters*rad ="<<4.605*emitx/pz<<" pi*meters*rad"<<std::endl;
  logger<<"*    90%emity="<< 4.605*mconstants::pi*emity/pz<<" meters*rad ="<<4.605*emity/pz<<" pi*meters*rad"<<std::endl;
  logger<<"*    90%emitz="<< 4.605*mconstants::pi*emitz*1.e9/(pconstants::c)<<" eV*s" <<std::endl;
  logger<<std::endl;
  logger<<std::endl;
  logger<<"*    95%emitx="<< 5.991*mconstants::pi*emitx/pz<<"  meters*rad ="<<5.991*emitx/pz<<" pi*meters*rad" <<std::endl;
  logger<<"*    95%emity="<< 5.991*mconstants::pi*emity/pz<<" meters*rad ="<<5.991*emity/pz<<" pi*meters*rad" <<std::endl;
  logger<<"*    95%emitz="<< 5.991*mconstants::pi*emitz*1.e9/(pconstants::c)<<" eV*s"<<std::endl;
  logger<<std::endl;
  logger<<"*    Normalized emitx="<< emitx*gamma*beta/pz<<" meters*rad ="<<emitx*gamma*beta/pz/mconstants::pi<<" pi*meters*rad"<<std::endl;
  logger<<"*    Normalized emity="<< emity*gamma*beta/pz<<" meters*rad ="<<emity*gamma*beta/pz/mconstants::pi<<" pi*meters*rad"<<std::endl;
  logger<<std::endl;
  logger<<"*    Normalized 90%emitx="<< 4.605*mconstants::pi*emitx*gamma*beta/pz<<"  meters*rad ="<<4.605*emitx*gamma*beta/pz<<" pi*meters*rad"<<std::endl;
  logger<<"*    Normalized 90%emity="<< 4.605*mconstants::pi*emity*gamma*beta/pz<<" meters*rad ="<<4.605*emity*gamma*beta/pz<<" pi*meters*rad"<<std::endl;
  logger<<std::endl;
  logger<<"*    Normalized 95%emitx="<< 5.991*mconstants::pi*emitx*gamma*beta/pz<<"  meters*rad ="<<5.991*emitx*gamma*beta/pz<<" pi*meters*rad"<<std::endl;
  logger<<"*    Normalized 95%emity="<< 5.991*mconstants::pi*emity*gamma*beta/pz<<" meters*rad ="<<5.991*emity*gamma*beta/pz<<" pi*meters*rad"<<std::endl;
  logger<<std::endl;
  logger<<"*    xrms="<<sqrt(mom2[0][0])/units[0] <<" meters"<<std::endl;
  logger<<"*    yrms="<<sqrt(mom2[2][2])/units[2]<<" meters"<<std::endl;
  logger<<"*    zrms="<<sqrt(mom2[4][4])/units[4] <<" meters="<<1e9*sqrt(mom2[4][4])/units[4]/pconstants::c/beta<<" ns  "<<std::endl;
  
  logger<<"*    pxrms="<<sqrt(mom2[1][1])/units[1] <<" GeV/c,  dpx/p="<<sqrt(mom2[1][1])<<std::endl;
  logger<<"*    pyrms="<<sqrt(mom2[3][3])/units[3] <<" GeV/c,   dpy/p="<<sqrt(mom2[3][3])<<std::endl;
  logger<<"*    pzrms="<<sqrt(mom2[5][5])/units[5] <<" GeV/c,  dpz/p="<<sqrt(mom2[5][5])<<std::endl;
  logger<<"*    Erms="<<sqrt(mom2[5][5])*beta*beta*energy<<" GeV,  deoe="<<sqrt(mom2[5][5])*beta*beta<<std::endl;
  
  logger<<"*    pz="<<pz<<"  GeV/c"<<std::endl;
  logger<<"*    total energy="<<energy<<"GeV,  kinetic energy="<<energy-pconstants::mp<<"GeV"<<std::endl;
  logger<<"****************************************************"<<std::endl;
#endif

} 
