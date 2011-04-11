#include "diagnostics.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "synergia/utils/eigen2/Eigen/Core"
#include "synergia/utils/eigen2/Eigen/LU"
#include <stdexcept>

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

MArray1d
Diagnostics::calculate_mean(Bunch const& bunch)
{
    MArray1d mean(boost::extents[6]);
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
    return mean;
}

MArray1d
Diagnostics::calculate_std(Bunch const& bunch, MArray1d_ref const& mean)
{
    MArray1d std(boost::extents[6]);
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
    return std;
}

MArray2d
Diagnostics::calculate_mom2(Bunch const& bunch, MArray1d_ref const& mean)
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

Diagnostics_basic::Diagnostics_basic(Bunch_sptr bunch_sptr,
        std::string const& filename) :
    bunch_sptr(bunch_sptr), filename(filename), have_writers(false), mean(
            boost::extents[6]), std(boost::extents[6])
{
}

bool
Diagnostics_basic::is_serial() const
{
    return true;
}

void
Diagnostics_basic::update()
{
    s = bunch_sptr->get_reference_particle().get_s();
    repetition = bunch_sptr->get_reference_particle().get_repetition();
    trajectory_length
            = bunch_sptr->get_reference_particle().get_trajectory_length();
    num_particles = bunch_sptr->get_total_num();
    real_num_particles = bunch_sptr->get_real_num();
    mean = Diagnostics::calculate_mean(*bunch_sptr);
    std = Diagnostics::calculate_std(*bunch_sptr, mean);
}

double
Diagnostics_basic::get_s() const
{
    return s;
}

int
Diagnostics_basic::get_repetition() const
{
    return repetition;
}

double
Diagnostics_basic::get_trajectory_length() const
{
    return trajectory_length;
}

int
Diagnostics_basic::get_num_particles() const
{
    return num_particles;
}

double
Diagnostics_basic::get_real_num_particles() const
{
    return real_num_particles;
}

Const_MArray1d_ref
Diagnostics_basic::get_mean() const
{
    return mean;
}

Const_MArray1d_ref
Diagnostics_basic::get_std() const
{
    return std;
}

void
Diagnostics_basic::init_writers(hid_t & hdf5_file)
{
    writer_s = new Hdf5_serial_writer<double > (hdf5_file, "s");
    writer_repetition = new Hdf5_serial_writer<int > (hdf5_file, "repetition");
    writer_trajectory_length = new Hdf5_serial_writer<double > (hdf5_file,
            "trajectory_length");
    writer_num_particles = new Hdf5_serial_writer<int > (hdf5_file,
            "num_particles");
    writer_real_num_particles = new Hdf5_serial_writer<double > (hdf5_file,
            "real_num_particles");
    writer_mean = new Hdf5_serial_writer<MArray1d_ref > (hdf5_file, "mean");
    writer_std = new Hdf5_serial_writer<MArray1d_ref > (hdf5_file, "std");
    have_writers = true;
}

void
Diagnostics_basic::write_hdf5()
{
    if (!have_writers) {
        throw(std::runtime_error(
                "Diagnostics_basic::write_hdf5 called before Diagnostics_basic::init_writers"));
    }
    writer_s->append(s);
    writer_repetition->append(repetition);
    writer_trajectory_length->append(trajectory_length);
    writer_num_particles->append(num_particles);
    writer_real_num_particles->append(real_num_particles);
    writer_mean->append(mean);
    writer_std->append(std);
}

Diagnostics_basic::~Diagnostics_basic()
{
    if (have_writers) {
        delete writer_s;
        delete writer_repetition;
        delete writer_trajectory_length;
        delete writer_num_particles;
        delete writer_real_num_particles;
        delete writer_mean;
        delete writer_std;
    }
}

void
Diagnostics_full2::update_full2()
{
    MArray2d sum2(boost::extents[6][6]);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            sum2[i][j] = 0.0;
        }
    }
    Const_MArray2d_ref particles(bunch_sptr->get_local_particles());
    for (int part = 0; part < bunch_sptr->get_local_num(); ++part) {
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
            bunch_sptr->get_comm().get());
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            mom2[i][j] = mom2[j][i] = mom2[i][j] / bunch_sptr->get_total_num();
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

Diagnostics_full2::Diagnostics_full2(Bunch_sptr bunch_sptr,
        std::string const& filename) :
    bunch_sptr(bunch_sptr), filename(filename), Diagnostics_basic(bunch_sptr,
            filename), have_writers(false), mom2(boost::extents[6][6]), corr(
            boost::extents[6][6])
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
    s = bunch_sptr->get_reference_particle().get_s();
    repetition = bunch_sptr->get_reference_particle().get_repetition();
    trajectory_length
            = bunch_sptr->get_reference_particle().get_trajectory_length();
    num_particles = bunch_sptr->get_total_num();
    real_num_particles = bunch_sptr->get_real_num();
    mean = calculate_mean(*bunch_sptr);
    update_full2();
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
    Diagnostics_basic::init_writers(hdf5_file);
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
    Diagnostics_basic::write_hdf5();
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

Diagnostics_particles::Diagnostics_particles(Bunch_sptr bunch_sptr,
        std::string const& filename, int max_particles) :
    bunch_sptr(bunch_sptr), filename(filename), max_particles(max_particles),
            have_writers(false)
{
}

bool
Diagnostics_particles::is_serial() const
{
    return false;
}

void
Diagnostics_particles::update()
{
}

void
Diagnostics_particles::init_writers(hid_t & hdf5_file)
{
    this->hdf5_file = hdf5_file;
    have_writers = true;
}

// jfa: this method is not complete! It doesn't work on multiple processors
void
Diagnostics_particles::write_hdf5()
{
    if (!have_writers) {
        throw(std::runtime_error(
                "Diagnostics_particles::write_hdf5 called before Diagnostics_particles::init_writers"));
    }
    Hdf5_writer<double > writer_pz(hdf5_file, "pz");
    double pz = bunch_sptr->get_reference_particle().get_momentum();
    writer_pz.write(pz);
    Hdf5_writer<double > writer_tlen(hdf5_file, "tlen");
    double tlen = bunch_sptr->get_reference_particle().get_trajectory_length();
    writer_tlen.write(tlen);
    Hdf5_writer<int > writer_rep(hdf5_file, "rep");
    int rep = bunch_sptr->get_reference_particle().get_repetition();
    writer_rep.write(rep);
    Hdf5_writer<double > writer_s(hdf5_file, "s");
    double s = bunch_sptr->get_reference_particle().get_s();
    writer_s.write(s);
    int local_num = bunch_sptr->get_local_num();
    Hdf5_chunked_array2d_writer
            writer_particles(hdf5_file, "particles",
                    bunch_sptr->get_local_particles()[boost::indices[range(0,
                            local_num)][range()]]);
    writer_particles.write_chunk(
            bunch_sptr->get_local_particles()[boost::indices[range(0, local_num)][range()]]);
}

Diagnostics_particles::~Diagnostics_particles()
{
}

Diagnostics_track::Diagnostics_track(Bunch_sptr bunch_sptr,
        std::string const& filename, int particle_id) :
    bunch_sptr(bunch_sptr), filename(filename), have_writers(false), coords(
            boost::extents[6]), found(false), particle_id(particle_id),
            last_index(-1)
{
}

bool
Diagnostics_track::is_serial() const
{
    return true;
}

void
Diagnostics_track::update()
{
    repetition = bunch_sptr->get_reference_particle().get_repetition();
    trajectory_length
            = bunch_sptr->get_reference_particle().get_trajectory_length();
    int index;
    found = false;
    if ((last_index > -1) && (last_index < bunch_sptr->get_local_num())) {
        if (particle_id
                == static_cast<int > (bunch_sptr->get_local_particles()[Bunch::id][last_index])) {
            index = last_index;
            found = true;
        }
    }
    if (!found) {
        index = 0;
        while ((particle_id
                != static_cast<int > (bunch_sptr->get_local_particles()[index][Bunch::id]))
                && (index < bunch_sptr->get_local_num())) {
            index += 1;
        }
        if (index < bunch_sptr->get_local_num()) {
            found = true;
        } else {
            found = false;
        }
    }
    if (found) {
        coords[0] = bunch_sptr->get_local_particles()[index][0];
        coords[1] = bunch_sptr->get_local_particles()[index][1];
        coords[2] = bunch_sptr->get_local_particles()[index][2];
        coords[3] = bunch_sptr->get_local_particles()[index][3];
        coords[4] = bunch_sptr->get_local_particles()[index][4];
        coords[5] = bunch_sptr->get_local_particles()[index][5];
        s = bunch_sptr->get_reference_particle().get_s();
        repetition = bunch_sptr->get_reference_particle().get_repetition();
        trajectory_length
                = bunch_sptr->get_reference_particle().get_trajectory_length();
    }
}

void
Diagnostics_track::init_writers(hid_t & hdf5_file)
{
    writer_coords = new Hdf5_serial_writer<MArray1d_ref > (hdf5_file, "coords");
    writer_s = new Hdf5_serial_writer<double > (hdf5_file, "s");
    writer_repetition = new Hdf5_serial_writer<int > (hdf5_file, "repetition");
    writer_trajectory_length = new Hdf5_serial_writer<double > (hdf5_file,
            "trajectory_length");
    have_writers = true;
}

void
Diagnostics_track::write_hdf5()
{
    if (!have_writers) {
        throw(std::runtime_error(
                "Diagnostics_track::write_hdf5 called before Diagnostics::init_writers"));
    }
    writer_coords->append(coords);
    writer_s->append(s);
    writer_repetition->append(repetition);
    writer_trajectory_length->append(trajectory_length);
}

Diagnostics_track::~Diagnostics_track()
{
    if (have_writers) {
        delete writer_coords;
        delete writer_s;
        delete writer_repetition;
        delete writer_trajectory_length;
    }
}
