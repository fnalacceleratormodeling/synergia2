#include "core_diagnostics.h"
#include <cmath>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

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
            sum[1] += lsum[0];
            sum[2] += lsum[0];
            sum[3] += lsum[0];
            sum[4] += lsum[0];
            sum[5] += lsum[0];
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
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 3; ++i) {
            sum[i] += particles[part][i * 2];
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
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 3; ++i) {
            double diff = particles[part][i * 2] - mean[i];
            sum[i] += diff * diff;
        }
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
