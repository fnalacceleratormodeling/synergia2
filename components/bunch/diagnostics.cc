#include "diagnostics.h"
#include <cmath>
#include "utils/eigen2/Eigen/Core"
#include "utils/eigen2/Eigen/LU"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

void
Diagnostics::update_mean(Bunch const& bunch)
{
    double sum[6] = { 0, 0, 0, 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            sum[i] += particles[part][i];
        }
    }
    MPI_Allreduce(sum, mean.origin(), 6, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 6; ++i) {
        mean[i] /= bunch.get_total_num();
    }
}

void
Diagnostics::update_std(Bunch const& bunch)
{
    double sum[6] = { 0, 0, 0, 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            double diff = particles[part][i] - mean[i];
            sum[i] += diff * diff;
        }
    }
    MPI_Allreduce(sum, std.origin(), 6, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 6; ++i) {
        std[i] = std::sqrt(std[i] / bunch.get_total_num());
    }
}

Diagnostics::Diagnostics() :
    mean(boost::extents[6]), std(boost::extents[6])
{
}

Diagnostics::Diagnostics(Bunch const& bunch, double s) :
    mean(boost::extents[6]), std(boost::extents[6])
{
    update(bunch, s);
}

void
Diagnostics::update(Bunch const& bunch, double s)
{
    this->s = s;
    update_mean(bunch);
    update_std(bunch);
}

double
Diagnostics::get_s() const
{
    return s;
}

Const_MArray1d_ref
Diagnostics::get_mean() const
{
    return mean;
}

Const_MArray1d_ref
Diagnostics::get_std() const
{
    return std;
}

Diagnostics::~Diagnostics()
{
}

void
Diagnostics_full2::update_full2(Bunch const& bunch)
{
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
        std[i] = std::sqrt(mom2[i][i]);
    }
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            corr[i][j] = corr[j][i] = mom2[i][j] / std::sqrt(mom2[i][i]
                    * mom2[j][j]);
        }
    }
}

void
Diagnostics_full2::update_emittances()
{
    Matrix<double, 6, 6 > mom2_matrix(mom2.origin());
    emitx = mom2_matrix.block<2, 2 > (Bunch::x, Bunch::x).determinant();
    emity = mom2_matrix.block<2, 2 > (Bunch::y, Bunch::y).determinant();
    emitz = mom2_matrix.block<2, 2 > (Bunch::z, Bunch::z).determinant();
    emitxy = mom2_matrix.block<4, 4 > (Bunch::x, Bunch::x).determinant();
    emitxyz = mom2_matrix.determinant();
}

Diagnostics_full2::Diagnostics_full2() :
    Diagnostics(), mom2(boost::extents[6][6]), corr(boost::extents[6][6])
{
}

Diagnostics_full2::Diagnostics_full2(Bunch const& bunch, double s) :
    Diagnostics(), mom2(boost::extents[6][6]), corr(boost::extents[6][6])
{
    update(bunch, s);
}

void
Diagnostics_full2::update(Bunch const& bunch, double s)
{
    this->s = s;
    update_mean(bunch);
    update_full2(bunch);
    update_emittances();
}

Const_MArray2d_ref
Diagnostics_full2::get_mom2() const
{
    return mom2;
}

Const_MArray2d_ref
Diagnostics_full2::get_corr() const
{
    return corr;
}

double
Diagnostics_full2::get_emitx() const
{
    return emitx;
}

double
Diagnostics_full2::get_emity() const
{
    return emity;
}

double
Diagnostics_full2::get_emitz() const
{
    return emitz;
}

double
Diagnostics_full2::get_emitxy() const
{
    return emitxy;
}

double
Diagnostics_full2::get_emitxyz() const
{
    return emitxyz;
}

Diagnostics_full2::~Diagnostics_full2()
{
}
