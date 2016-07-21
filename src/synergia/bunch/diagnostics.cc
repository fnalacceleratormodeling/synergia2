#include "diagnostics.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "Eigen/Core"
#include "Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

using namespace Eigen;

Diagnostics ::Diagnostics(std::string const& name, std::string const& filename,
        std::string const& local_dir) :
        name(name), filename(filename), local_dir(local_dir), have_bunch_(
                false), write_helper_ptr(0), have_write_helper_(0)
{
}

std::string const&
Diagnostics::get_filename() const
{
    return filename;
}

std::string const&
Diagnostics::get_local_dir() const
{
    return local_dir;
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
            is_serial(), get_bunch().get_comm_sptr(), local_dir);
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

template<class Archive>
    void
    Diagnostics::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(name);
        ar & BOOST_SERIALIZATION_NVP(filename);
        ar & BOOST_SERIALIZATION_NVP(local_dir);
        ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
        ar & BOOST_SERIALIZATION_NVP(have_bunch_);
        ar & BOOST_SERIALIZATION_NVP(write_helper_ptr);
        ar & BOOST_SERIALIZATION_NVP(have_write_helper_);
    }

template
void
Diagnostics::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics::~Diagnostics()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics)
