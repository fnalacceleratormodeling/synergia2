#include "diagnostics_full2.h"
#include <cmath>
#include "synergia/utils/eigen2/Eigen/Core"
#include "synergia/utils/eigen2/Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

const char Diagnostics_full2::name[] = "diagnostics_full2";

void
Diagnostics_full2::update_full2()
{
    MArray2d sum2(boost::extents[6][6]);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            sum2[i][j] = 0.0;
        }
    }
    Const_MArray2d_ref particles(get_bunch().get_local_particles());
    for (int part = 0; part < get_bunch().get_local_num(); ++part) {
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
            get_bunch().get_comm().get());
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            mom2[i][j] = mom2[j][i] = mom2[i][j] / get_bunch().get_total_num();
        }
        std[i] = std::sqrt(mom2[i][i]);
    }
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            corr[i][j] = corr[j][i] = mom2[i][j] / std::sqrt(
                    mom2[i][i] * mom2[j][j]);
        }
    }
}

void
Diagnostics_full2::update_emittances()
{
    Matrix<double, 6, 6 > mom2_matrix(mom2.origin());
    emitx = std::sqrt(
            mom2_matrix.block<2, 2 > (Bunch::x, Bunch::x).determinant());
    emity = std::sqrt(
            mom2_matrix.block<2, 2 > (Bunch::y, Bunch::y).determinant());
    emitz = std::sqrt(
            mom2_matrix.block<2, 2 > (Bunch::z, Bunch::z).determinant());
    emitxy = std::sqrt(
            mom2_matrix.block<4, 4 > (Bunch::x, Bunch::x).determinant());
    emitxyz = std::sqrt(mom2_matrix.determinant());
}

Diagnostics_full2::Diagnostics_full2(std::string const& filename) :
    Diagnostics_full2::Diagnostics(Diagnostics_full2::name, filename),
            have_writers(false), mean(boost::extents[6]),
            std(boost::extents[6]),
            min(boost::extents[3]), max(boost::extents[3]),
            mom2(boost::extents[6][6]),
            corr(boost::extents[6][6]), writer_s(0), writer_repetition(0),
            writer_trajectory_length(0), writer_num_particles(0),
            writer_real_num_particles(0), writer_mean(0), writer_std(0),
            writer_min(0), writer_max(0),
            writer_mom2(0), writer_corr(0), writer_emitx(0), writer_emity(0),
            writer_emitz(0), writer_emitxy(0), writer_emitxyz(0)

{
}

Diagnostics_full2::Diagnostics_full2()
{
}

bool
Diagnostics_full2::is_serial() const
{
    return true;
}

void
Diagnostics_full2::update()
{
    get_bunch().convert_to_state(Bunch::fixed_z_lab);
    s = get_bunch().get_reference_particle().get_s();
    repetition = get_bunch().get_reference_particle().get_repetition();
    trajectory_length
            = get_bunch().get_reference_particle().get_trajectory_length();
    num_particles = get_bunch().get_total_num();
    real_num_particles = get_bunch().get_real_num();
    min = Core_diagnostics::calculate_min(get_bunch());
    max = Core_diagnostics::calculate_max(get_bunch());
    mean = Core_diagnostics::calculate_mean(get_bunch());
    update_full2();
    update_emittances();
}

double
Diagnostics_full2::get_s() const
{
    return s;
}

int
Diagnostics_full2::get_repetition() const
{
    return repetition;
}

double
Diagnostics_full2::get_trajectory_length() const
{
    return trajectory_length;
}

int
Diagnostics_full2::get_num_particles() const
{
    return num_particles;
}

double
Diagnostics_full2::get_real_num_particles() const
{
    return real_num_particles;
}

Const_MArray1d_ref
Diagnostics_full2::get_mean() const
{
    return mean;
}

Const_MArray1d_ref
Diagnostics_full2::get_std() const
{
    return std;
}

const MArray1d
Diagnostics_full2::get_min() const
{
    return min;
}

const MArray1d
Diagnostics_full2::get_max() const
{
    return max;
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
Diagnostics_full2::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");
        writer_repetition = new Hdf5_serial_writer<int > (file_sptr,
                "repetition");
        writer_trajectory_length = new Hdf5_serial_writer<double > (file_sptr,
                "trajectory_length");
        writer_num_particles = new Hdf5_serial_writer<int > (file_sptr,
                "num_particles");
        writer_real_num_particles = new Hdf5_serial_writer<double > (file_sptr,
                "real_num_particles");
        writer_mean = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "mean");
        writer_std = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "std");
        writer_min = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "min");
        writer_max = new Hdf5_serial_writer<MArray1d_ref > (file_sptr, "max");
        writer_mom2 = new Hdf5_serial_writer<MArray2d_ref > (file_sptr, "mom2");
        writer_corr = new Hdf5_serial_writer<MArray2d_ref > (file_sptr, "corr");
        writer_emitx = new Hdf5_serial_writer<double > (file_sptr, "emitx");
        writer_emity = new Hdf5_serial_writer<double > (file_sptr, "emity");
        writer_emitz = new Hdf5_serial_writer<double > (file_sptr, "emitz");
        writer_emitxy = new Hdf5_serial_writer<double > (file_sptr, "emitxy");
        writer_emitxyz = new Hdf5_serial_writer<double > (file_sptr, "emitxyz");
        have_writers = true;
    }
}

void
Diagnostics_full2::write()
{
    if (get_write_helper().write_locally()) {
        init_writers(get_write_helper().get_hdf5_file_sptr());
        writer_s->append(s);
        writer_repetition->append(repetition);
        writer_trajectory_length->append(trajectory_length);
        writer_num_particles->append(num_particles);
        writer_real_num_particles->append(real_num_particles);
        writer_mean->append(mean);
        writer_std->append(std);
        writer_min->append(min);
        writer_max->append(max);
        writer_mom2->append(mom2);
        writer_corr->append(corr);
        writer_emitx->append(emitx);
        writer_emity->append(emity);
        writer_emitz->append(emitz);
        writer_emitxy->append(emitxy);
        writer_emitxyz->append(emitxyz);
        get_write_helper().finish_write();
    }
}

Diagnostics_full2::~Diagnostics_full2()
{
    if (have_writers) {
        delete writer_emitxyz;
        delete writer_emitxy;
        delete writer_emitz;
        delete writer_emity;
        delete writer_emitx;
        delete writer_corr;
        delete writer_mom2;
        delete writer_max;
        delete writer_min;
        delete writer_std;
        delete writer_mean;
        delete writer_real_num_particles;
        delete writer_num_particles;
        delete writer_trajectory_length;
        delete writer_repetition;
        delete writer_s;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_full2)
