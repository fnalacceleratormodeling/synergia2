#include "diagnostics_basic.h"
#include "core_diagnostics.h"
#include <cmath>

const char Diagnostics_basic::name[] = "diagnostics_basic";

Diagnostics_basic::Diagnostics_basic(std::string const& filename) :
    Diagnostics_basic::Diagnostics(Diagnostics_basic::name, filename),
            have_writers(false), mean(boost::extents[6]),
            std(boost::extents[6]),
            min(boost::extents[3]), max(boost::extents[3]),
            writer_s(0), writer_repetition(0),
            writer_trajectory_length(0), writer_num_particles(0),
            writer_real_num_particles(0), writer_mean(0), writer_std(0),
            writer_min(0), writer_max(0)
{
}

Diagnostics_basic::Diagnostics_basic()
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
    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
    s = get_bunch().get_reference_particle().get_s();
    repetition = get_bunch().get_reference_particle().get_repetition();
    trajectory_length
            = get_bunch().get_reference_particle().get_trajectory_length();
    num_particles = get_bunch().get_total_num();
    real_num_particles = get_bunch().get_real_num();
    mean = Core_diagnostics::calculate_mean(get_bunch());
    std = Core_diagnostics::calculate_std(get_bunch(), mean);
    min = Core_diagnostics::calculate_min(get_bunch());
    max = Core_diagnostics::calculate_max(get_bunch());
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

const MArray1d
Diagnostics_basic::get_min() const
{
    return min;
}

const MArray1d
Diagnostics_basic::get_max() const
{
    return max;
}

void
Diagnostics_basic::init_writers(Hdf5_file_sptr file_sptr)
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
        have_writers = true;
    }
}

void
Diagnostics_basic::write()
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
        get_write_helper().finish_write();
    }
}

Diagnostics_basic::~Diagnostics_basic()
{
    if (have_writers) {
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
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_basic)
