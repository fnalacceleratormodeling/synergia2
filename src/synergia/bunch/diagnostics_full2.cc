#include "diagnostics_full2.h"
#include <cmath>
#include "Eigen/Core"
#include "Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

using namespace Eigen;

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

namespace
{
  const double tiny = 1.0e-15;
  double
  eliminate_small_negative(double x)
  {
    if ((x < 0) && (x > -tiny)) {
      return 0;
    } else {
      return x;
    }
  }
}

void
Diagnostics_full2::update_emittances()
{
    Matrix<double, 6, 6 > mom2_matrix(mom2.origin());
    emitx = std::sqrt(eliminate_small_negative(
            mom2_matrix.block<2, 2 > (Bunch::x, Bunch::x).determinant()));
    emity = std::sqrt(eliminate_small_negative(
            mom2_matrix.block<2, 2 > (Bunch::y, Bunch::y).determinant()));
    emitz = std::sqrt(eliminate_small_negative(
            mom2_matrix.block<2, 2 > (Bunch::z, Bunch::z).determinant()));
    emitxy = std::sqrt(eliminate_small_negative(
            mom2_matrix.block<4, 4 > (Bunch::x, Bunch::x).determinant()));
    emitxyz = std::sqrt(eliminate_small_negative(mom2_matrix.determinant()));
}

Diagnostics_full2::Diagnostics_full2(std::string const& filename,
        std::string const& local_dir) :
        Diagnostics_full2::Diagnostics(Diagnostics_full2::name, filename,
                local_dir), have_writers(false), writer_s_n(0), writer_repetition(
                0), writer_s(0), writer_num_particles(0), writer_real_num_particles(
                0), writer_pz(0), mean(boost::extents[6]), writer_mean(0), std(
                boost::extents[6]), writer_std(0), min(boost::extents[3]), writer_min(
                0), max(boost::extents[3]), writer_max(0), mom2(
                boost::extents[6][6]), writer_mom2(0), corr(
                boost::extents[6][6]), writer_corr(0), writer_emitx(0), writer_emity(
                0), writer_emitz(0), writer_emitxy(0), writer_emitxyz(0)

{
}

Diagnostics_full2::Diagnostics_full2() : have_writers(false)
{
}

bool
Diagnostics_full2::is_serial() const
{
    return true;
}

void
Diagnostics_full2::update()
{   if (get_bunch().get_comm().has_this_rank()){
      get_bunch().convert_to_state(Bunch::fixed_z_lab);
      s_n = get_bunch().get_reference_particle().get_s_n();
      repetition = get_bunch().get_reference_particle().get_repetition();
      s
	      = get_bunch().get_reference_particle().get_s();
      num_particles = get_bunch().get_total_num();
      real_num_particles = get_bunch().get_real_num();
      pz = get_bunch().get_reference_particle().get_momentum();
      min = Core_diagnostics::calculate_min(get_bunch());
      max = Core_diagnostics::calculate_max(get_bunch());
      mean = Core_diagnostics::calculate_mean(get_bunch());
      update_full2();
      update_emittances();
    }
}

double
Diagnostics_full2::get_s_n() const
{
    return s_n;
}

int
Diagnostics_full2::get_repetition() const
{
    return repetition;
}

double
Diagnostics_full2::get_s() const
{
    return s;
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

double
Diagnostics_full2::get_pz() const
{
    return pz;
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
        Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );
        int chg = get_bunch().get_reference_particle().get_charge();
        file_sptr->write(chg, "charge");
        double pmass = fourp.get_mass();
        file_sptr->write(pmass, "mass");
        writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
        writer_repetition = new Hdf5_serial_writer<int > (file_sptr,
                "repetition");
        writer_s = new Hdf5_serial_writer<double > (file_sptr,
                "s");
        writer_num_particles = new Hdf5_serial_writer<int > (file_sptr,
                "num_particles");
        writer_real_num_particles = new Hdf5_serial_writer<double > (file_sptr,
                "real_num_particles");
        writer_pz = new Hdf5_serial_writer<double > (file_sptr,"pz");
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
    if (get_bunch().get_comm().has_this_rank()){
      if (get_write_helper().write_locally()) {
	  init_writers(get_write_helper().get_hdf5_file_sptr());
	  writer_s_n->append(s_n);
	  writer_repetition->append(repetition);
	  writer_s->append(s);
	  writer_num_particles->append(num_particles);
	  writer_real_num_particles->append(real_num_particles);
      writer_pz->append(pz);
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
}

// begin egs screwing around
bool Diagnostics_full2::get_have_writers() { return have_writers; }
// end egs screwing around
template<class Archive>
    void
    Diagnostics_full2::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                & BOOST_SERIALIZATION_NVP(have_writers)
                & BOOST_SERIALIZATION_NVP(s_n)
                & BOOST_SERIALIZATION_NVP(writer_s_n)
                & BOOST_SERIALIZATION_NVP(repetition)
                & BOOST_SERIALIZATION_NVP(writer_repetition)
                & BOOST_SERIALIZATION_NVP(s)
                & BOOST_SERIALIZATION_NVP(writer_s)
                & BOOST_SERIALIZATION_NVP(num_particles)
                & BOOST_SERIALIZATION_NVP(writer_num_particles)
                & BOOST_SERIALIZATION_NVP(real_num_particles)
                & BOOST_SERIALIZATION_NVP(writer_real_num_particles)
                & BOOST_SERIALIZATION_NVP(pz)
                & BOOST_SERIALIZATION_NVP(writer_pz)
                & BOOST_SERIALIZATION_NVP(mean)
                & BOOST_SERIALIZATION_NVP(writer_mean)
                & BOOST_SERIALIZATION_NVP(std)
                & BOOST_SERIALIZATION_NVP(writer_std)
                & BOOST_SERIALIZATION_NVP(min)
                & BOOST_SERIALIZATION_NVP(writer_min)
                & BOOST_SERIALIZATION_NVP(max)
                & BOOST_SERIALIZATION_NVP(writer_max)
                & BOOST_SERIALIZATION_NVP(mom2)
                & BOOST_SERIALIZATION_NVP(writer_mom2)
                & BOOST_SERIALIZATION_NVP(corr)
                & BOOST_SERIALIZATION_NVP(writer_corr)
                & BOOST_SERIALIZATION_NVP(emitx)
                & BOOST_SERIALIZATION_NVP(emity)
                & BOOST_SERIALIZATION_NVP(emitz)
                & BOOST_SERIALIZATION_NVP(emitxy)
                & BOOST_SERIALIZATION_NVP(emitxyz)
                & BOOST_SERIALIZATION_NVP(writer_emitx)
                & BOOST_SERIALIZATION_NVP(writer_emity)
                & BOOST_SERIALIZATION_NVP(writer_emitz)
                & BOOST_SERIALIZATION_NVP(writer_emitxy)
                & BOOST_SERIALIZATION_NVP(writer_emitxyz);
    }

template
void
Diagnostics_full2::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_full2::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_full2::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_full2::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

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
        delete writer_pz;
        delete writer_real_num_particles;
        delete writer_num_particles;
        delete writer_s;
        delete writer_repetition;
        delete writer_s_n;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_full2)
