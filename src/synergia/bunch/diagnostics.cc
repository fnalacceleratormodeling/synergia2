#include "diagnostics.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "synergia/utils/eigen2/Eigen/Core"
#include "synergia/utils/eigen2/Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

Diagnostics::Diagnostics(std::string const& name, std::string const& filename) :
    name(name), filename(filename), have_bunch_(false), write_helper_ptr(0),
            have_write_helper_(0)
{
}

std::string const&
Diagnostics::get_filename() const
{
    return filename;
}

void
Diagnostics::set_bunch_sptr(Bunch_sptr bunch_sptr)
{
    this->bunch_sptr = bunch_sptr;
    have_bunch_ = true;
}

bool
Diagnostics::have_bunch() const
{
    return have_bunch_;
}

Bunch &
Diagnostics::get_bunch()
{
    if (!have_bunch_) {
        throw std::runtime_error(name + ": bunch not set");
    }
    return *bunch_sptr;
}

void
Diagnostics::delete_write_helper_ptr()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
        have_write_helper_ = false;
    }
}

Diagnostics_write_helper *
Diagnostics::new_write_helper_ptr()
{
    delete_write_helper_ptr();
    return new Diagnostics_write_helper(get_filename(),
            is_serial(), get_bunch().get_comm());
}

bool
Diagnostics::have_write_helper() const
{
    return have_write_helper_;
}

Diagnostics_write_helper &
Diagnostics::get_write_helper()
{
    if (!have_write_helper_) {
        write_helper_ptr = new_write_helper_ptr();
        have_write_helper_ = true;
    }
    return *write_helper_ptr;
}

Diagnostics::Diagnostics()
{
}

Diagnostics::~Diagnostics()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
    }
}

BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics)




const char Diagnostics_reference_particle::name[] = "diagnostics_reference_particle";

Diagnostics_reference_particle::Diagnostics_reference_particle(
        std::string const& filename) :
            Diagnostics_reference_particle::Diagnostics(
                    Diagnostics_reference_particle::name, filename),
            have_writers(false), writer_beta(0), writer_gamma(0),
            writer_state(0), writer_s(0)
{
}

Diagnostics_reference_particle::Diagnostics_reference_particle()
{
}

bool
Diagnostics_reference_particle::is_serial() const
{
    return true;
}

void
Diagnostics_reference_particle::update()
{
}

void
Diagnostics_reference_particle::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_beta = new Hdf5_serial_writer<double > (file_sptr, "beta");
        writer_gamma = new Hdf5_serial_writer<double > (file_sptr, "gamma");
        writer_state = new Hdf5_serial_writer<MArray1d_ref > (file_sptr,
                "state");
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");
        have_writers = true;
    }
}

void
Diagnostics_reference_particle::write()
{
    if (get_write_helper().write_locally()) {
        init_writers(get_write_helper().get_hdf5_file_sptr());
        double beta = get_bunch().get_reference_particle().get_beta();
        writer_beta->append(beta);
        double gamma = get_bunch().get_reference_particle().get_gamma();
        writer_gamma->append(gamma);
        MArray1d state(get_bunch().get_reference_particle().get_state());
        writer_state->append(state);
        double s = get_bunch().get_reference_particle().get_s();
        writer_s->append(s);
        get_write_helper().finish_write();
    }
}

Diagnostics_reference_particle::~Diagnostics_reference_particle()
{
    if (have_writers) {
        delete writer_beta;
        delete writer_gamma;
        delete writer_state;
        delete writer_s;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_reference_particle)

