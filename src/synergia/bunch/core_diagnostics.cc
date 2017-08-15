#include "core_diagnostics.h"
#include <cmath>
#include "Eigen/Core"
#include "Eigen/LU"
#include <stdexcept>
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/synergia_omp.h"

using namespace Eigen;

MArray1d
Core_diagnostics::calculate_mean(Bunch const& bunch)
{
    MArray1d mean(boost::extents[6]);
    double sum[6] = { 0, 0, 0, 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    int npart = bunch.get_local_num();

    #pragma omp parallel shared(npart, particles)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int l = npart / nt;
        int s = it * l;
        int e = (it==nt-1) ? npart : (it+1)*l;

        double lsum[6] = { 0, 0, 0, 0, 0, 0 };

        for (int part = s; part < e; ++part) 
        {
            lsum[0] += particles[part][0];
            lsum[1] += particles[part][1];
            lsum[2] += particles[part][2];
            lsum[3] += particles[part][3];
            lsum[4] += particles[part][4];
            lsum[5] += particles[part][5];
        }

        #pragma omp critical
        {
            sum[0] += lsum[0];
            sum[1] += lsum[1];
            sum[2] += lsum[2];
            sum[3] += lsum[3];
            sum[4] += lsum[4];
            sum[5] += lsum[5];
        } // end of omp critical

    } // end of omp parallel

    double t;
    t = simple_timer_current();
    MPI_Allreduce(sum, mean.origin(), 6, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    t = simple_timer_show(t, "allmpireduce_in_diagnostic mean");

    for (int i = 0; i < 6; ++i) {
        mean[i] /= bunch.get_total_num();
    }
    return mean;
}

double
Core_diagnostics::calculate_z_mean(Bunch const& bunch)
{
    double sum = 0;
    double mean;
    Const_MArray2d_ref particles(bunch.get_local_particles());
    int npart = bunch.get_local_num();

    #pragma omp parallel shared(npart, particles)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int l = npart / nt;
        int s = it * l;
        int e = (it==nt-1) ? npart : (it+1)*l;

        double lsum = 0;

        for (int part = s; part < e; ++part) 
        {
            lsum += particles[part][4];
        }

        #pragma omp critical
        sum += lsum;

    } // end of omp parallel

    MPI_Allreduce(&sum, &mean, 1, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    mean /= bunch.get_total_num();
    return mean;
}

double
Core_diagnostics::calculate_z_std(Bunch const& bunch, double const& mean)
{
    double sum = 0;
    double std;
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double diff = particles[part][4] - mean;
        sum += diff * diff;
    }
    MPI_Allreduce(&sum, &std, 1, MPI_DOUBLE, MPI_SUM, bunch.get_comm().get());
    std = std::sqrt(std / bunch.get_total_num());
    return std;
}

MArray1d
Core_diagnostics::calculate_spatial_mean(Bunch const& bunch)
{
    MArray1d mean(boost::extents[3]);
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
    return mean;
}

MArray1d
Core_diagnostics::calculate_std(Bunch const& bunch, MArray1d_ref const& mean)
{
    MArray1d std(boost::extents[6]);
    double sum[6] = { 0, 0, 0, 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    int npart = bunch.get_local_num();

    #pragma omp parallel shared(npart, particles)
    {
        int nt = omp_get_num_threads();
        int it = omp_get_thread_num();

        int l = npart / nt;
        int s = it * l;
        int e = (it==nt-1) ? npart : (it+1)*l;

        double lsum[6] = { 0, 0, 0, 0, 0, 0 };
        double diff;

        for(int part = s; part < e; ++part)
        {
            diff = particles[part][0] - mean[0]; lsum[0] += diff * diff;
            diff = particles[part][1] - mean[1]; lsum[1] += diff * diff;
            diff = particles[part][2] - mean[2]; lsum[2] += diff * diff;
            diff = particles[part][3] - mean[3]; lsum[3] += diff * diff;
            diff = particles[part][4] - mean[4]; lsum[4] += diff * diff;
            diff = particles[part][5] - mean[5]; lsum[5] += diff * diff;
        }
 
        #pragma omp critical
        {
            sum[0] += lsum[0];
            sum[1] += lsum[1];
            sum[2] += lsum[2];
            sum[3] += lsum[3];
            sum[4] += lsum[4];
            sum[5] += lsum[5];
        } // end of omp critical
    }

    MPI_Allreduce(sum, std.origin(), 6, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());

    for (int i = 0; i < 6; ++i) {
        std[i] = std::sqrt(std[i] / bunch.get_total_num());
    }

    return std;
}

MArray1d
Core_diagnostics::calculate_spatial_std(Bunch const& bunch, MArray1d_ref const& mean)
{
    MArray1d std(boost::extents[3]);
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
    return std;
}

MArray2d
Core_diagnostics::calculate_mom2(Bunch const& bunch, MArray1d_ref const& mean)
{
    MArray2d mom2(boost::extents[6][6]);
    MArray2d sum2(boost::extents[6][6]);

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
    return mom2;
}

MArray1d
Core_diagnostics::calculate_min(Bunch const& bunch)
{
    MArray1d min(boost::extents[3]);
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

    return min;
}

MArray1d
Core_diagnostics::calculate_max(Bunch const& bunch)
{
    MArray1d max(boost::extents[3]);
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

    return max;
}

void
Core_diagnostics::print_bunch_parameters(MArray2d_ref const& mom2, double beta)
{
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

} 
