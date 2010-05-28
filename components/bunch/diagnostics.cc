#include "diagnostics.h"
#include "utils/hdf5_writer.h"
#include "utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "utils/eigen2/Eigen/Core"
#include "utils/eigen2/Eigen/LU"
#include <stdexcept>

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
    mean(boost::extents[6]), std(boost::extents[6]), have_writers(false)
{
}

Diagnostics::Diagnostics(Bunch const& bunch) :
    mean(boost::extents[6]), std(boost::extents[6]), have_writers(false)
{
    update(bunch);
}

void
Diagnostics::update(Bunch const& bunch)
{
    s = bunch.get_reference_particle().get_s();
    repetition = bunch.get_reference_particle().get_repetition();
    trajectory_length = bunch.get_reference_particle().get_trajectory_length();
    update_mean(bunch);
    update_std(bunch);
}

double
Diagnostics::get_s() const
{
    return s;
}

int
Diagnostics::get_repetition() const
{
    return repetition;
}

double
Diagnostics::get_trajectory_length() const
{
    return trajectory_length;
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

void
Diagnostics::init_writers(hid_t & hdf5_file)
{
    writer_s = new Hdf5_serial_writer<double > (hdf5_file, "s");
    writer_repetition = new Hdf5_serial_writer<int > (hdf5_file, "repetition");
    writer_trajectory_length = new Hdf5_serial_writer<double > (hdf5_file,
            "trajectory_length");
    writer_mean = new Hdf5_serial_writer<MArray1d_ref > (hdf5_file, "mean");
    writer_std = new Hdf5_serial_writer<MArray1d_ref > (hdf5_file, "std");
    have_writers = true;
}

void
Diagnostics::write_hdf5()
{
    if (!have_writers) {
        throw(std::runtime_error(
                "Diagnostics::write_hdf5 called before Diagnostics::init_writers"));
    }
    writer_s->append(s);
    writer_repetition->append(repetition);
    writer_trajectory_length->append(trajectory_length);
    writer_mean->append(mean);
    writer_std->append(std);
}

Diagnostics::~Diagnostics()
{
    if (have_writers) {
        delete writer_s;
        delete writer_repetition;
        delete writer_trajectory_length;
        delete writer_mean;
        delete writer_std;
    }
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
    Diagnostics(), mom2(boost::extents[6][6]), corr(boost::extents[6][6]),
            have_writers(false)
{
}

Diagnostics_full2::Diagnostics_full2(Bunch const& bunch) :
    Diagnostics(), mom2(boost::extents[6][6]), corr(boost::extents[6][6]),
            have_writers(false)
{
    update(bunch);
}

void
Diagnostics_full2::update(Bunch const& bunch)
{
    s = bunch.get_reference_particle().get_s();
    repetition = bunch.get_reference_particle().get_repetition();
    trajectory_length = bunch.get_reference_particle().get_trajectory_length();
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

void
Diagnostics_full2::init_writers(hid_t & hdf5_file)
{
    Diagnostics::init_writers(hdf5_file);
    writer_mom2 = new Hdf5_serial_writer<MArray2d_ref > (hdf5_file, "mom2");
    writer_corr = new Hdf5_serial_writer<MArray2d_ref > (hdf5_file, "corr");
    writer_emitx = new Hdf5_serial_writer<double > (hdf5_file, "emitx");
    writer_emity = new Hdf5_serial_writer<double > (hdf5_file, "emity");
    writer_emitz = new Hdf5_serial_writer<double > (hdf5_file, "emitz");
    writer_emitxy = new Hdf5_serial_writer<double > (hdf5_file, "emitxy");
    writer_emitxyz = new Hdf5_serial_writer<double > (hdf5_file, "emitxyz");
    have_writers = true;
}

void
Diagnostics_full2::write_hdf5()
{
    Diagnostics::write_hdf5();
    writer_mom2->append(mom2);
    writer_corr->append(corr);
    writer_emitx->append(emitx);
    writer_emity->append(emity);
    writer_emitz->append(emitz);
    writer_emitxy->append(emitxy);
    writer_emitxyz->append(emitxyz);
}

Diagnostics_full2::~Diagnostics_full2()
{
    if (have_writers) {
        delete writer_mom2;
        delete writer_corr;
        delete writer_emitx;
        delete writer_emity;
        delete writer_emitz;
        delete writer_emitxy;
        delete writer_emitxyz;
    }
}

Diagnostics_particles::Diagnostics_particles(int max_particles) :
    bunch_ptr(0), max_particles(max_particles)
{
}

Diagnostics_particles::Diagnostics_particles(Bunch const& bunch,
        int max_particles) :
    bunch_ptr(&bunch), max_particles(max_particles)
{
}

void
Diagnostics_particles::update(Bunch const& bunch)
{
    bunch_ptr = &bunch;
}

void
Diagnostics_particles::init_writers(hid_t & hdf5_file)
{
    this->hdf5_file = hdf5_file;
}

// jfa: this method is not complete! It doesn't work on multiple processors
// or implement the max_particles parameter.
void
Diagnostics_particles::write_hdf5()
{
    Hdf5_writer<double > writer_pz(hdf5_file, "pz");
    double pz = bunch_ptr->get_reference_particle().get_momentum();
    writer_pz.write(pz);
    Hdf5_chunked_array2d_writer writer_particles(hdf5_file, "particles",
            bunch_ptr->get_local_particles());
    writer_particles.write_chunk(bunch_ptr->get_local_particles());
}

Diagnostics_particles::~Diagnostics_particles()
{
}
