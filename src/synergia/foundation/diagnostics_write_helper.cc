#include "diagnostics_write_helper.h"
#include <sstream>
#include <iomanip>
#include <stdexcept>

// avoid bad interaction between Boost Filesystem and clang
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>

void
move_file_overwrite_if_exists(std::string const & source,
        std::string const & dest)
{
    if (boost::filesystem::exists(dest)) {
        boost::filesystem::remove(dest);
    }
    boost::filesystem::copy_file(source, dest);
    boost::filesystem::remove(source);
}

std::string
Diagnostics_write_helper::get_filename(bool include_local_dir)
{
    std::stringstream sstream;
    if (include_local_dir && (local_dir != "")) {
        sstream << local_dir;
        sstream << "/";
    }
    sstream << filename_base;
    if (!serial) {
        sstream << "_";
        sstream << std::setw(4);
        sstream << std::setfill('0');
        sstream << count;
    }
    sstream << filename_suffix;
    return sstream.str();
}

void
Diagnostics_write_helper::open_file()
{
    if (write_locally() && !have_file) {
        file_sptr = Hdf5_file_sptr(
                new Hdf5_file(get_filename(true).c_str(), Hdf5_file::truncate));
        have_file = true;
    }
}

Diagnostics_write_helper::Diagnostics_write_helper(std::string const& filename,
		bool serial, Commxx_sptr commxx_sptr, std::string const& local_dir,
		int writer_rank) :
		filename(filename), local_dir(local_dir), serial(serial), commxx_sptr(
				commxx_sptr), have_file(false), count(0) {
	if (writer_rank == default_rank) {
		this->writer_rank = commxx_sptr->get_size() - 1;
    } else {
        this->writer_rank = writer_rank;
    }
    std::string::size_type idx = filename.rfind('.');
    if (idx == std::string::npos) {
        filename_base = filename;
        filename_suffix = "";
    } else {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }
}

Diagnostics_write_helper::Diagnostics_write_helper()
{
}

int
Diagnostics_write_helper::get_count() const
{
    return count;
}

void
Diagnostics_write_helper::set_count(int count)
{
    this->count = count;
}

void
Diagnostics_write_helper::increment_count()
{
    ++count;
}

bool
Diagnostics_write_helper::write_locally()
{
    return commxx_sptr->get_rank() == writer_rank;
}

int
Diagnostics_write_helper::get_writer_rank()
{
    return writer_rank;
}

Hdf5_file_sptr
Diagnostics_write_helper::get_hdf5_file_sptr()
{
    if (!write_locally()) {
        throw std::runtime_error(
                "Diagnostics_write_helper::get_hdf5_file_sptr() called on a non-writer rank.");
    }
    open_file();
    return file_sptr;
}

void
Diagnostics_write_helper::finish_write()
{
    if (write_locally() && !serial) {
        file_sptr->close();
        if (local_dir != "") {
            move_file_overwrite_if_exists(get_filename(true), get_filename(false));
        }
        file_sptr.reset();
        have_file = false;
    }
    ++count;
    if (write_locally() && serial) {
        if (count % flush_period == 0) {
            file_sptr->flush();
        }
    }
}

template<class Archive>
    void
    Diagnostics_write_helper::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(writer_rank)
                & BOOST_SERIALIZATION_NVP(filename)
                & BOOST_SERIALIZATION_NVP(local_dir)
                & BOOST_SERIALIZATION_NVP(serial)
                & BOOST_SERIALIZATION_NVP(commxx_sptr)
                & BOOST_SERIALIZATION_NVP(file_sptr)
                & BOOST_SERIALIZATION_NVP(have_file)
                & BOOST_SERIALIZATION_NVP(count)
                & BOOST_SERIALIZATION_NVP(filename_base)
                & BOOST_SERIALIZATION_NVP(filename_suffix);
    }

template
void
Diagnostics_write_helper::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_write_helper::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_write_helper::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_write_helper::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_write_helper::~Diagnostics_write_helper()
{
    if (write_locally() && serial && have_file) {
        file_sptr->close();
        if (local_dir != "") {
            move_file_overwrite_if_exists(get_filename(true),
                    get_filename(false));
        }
    }
}
