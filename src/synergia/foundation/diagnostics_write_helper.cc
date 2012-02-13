#include "diagnostics_write_helper.h"
#include <sstream>
#include <iomanip>
#include <stdexcept>

std::string
Diagnostics_write_helper::get_filename()
{
    std::stringstream sstream;
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
        file_sptr = boost::shared_ptr<H5::H5File >(
                new H5::H5File(get_filename().c_str(), H5F_ACC_TRUNC));
        have_file = true;
    }
}

void
Diagnostics_write_helper::reopen_file()
{
    if (write_locally() && !have_file) {
        file_sptr = boost::shared_ptr<H5::H5File >(
                new H5::H5File(get_filename().c_str(), H5F_ACC_RDWR));
        have_file = true;
    }
}

void
Diagnostics_write_helper::construct(std::string const& filename, bool serial,
        int write_skip, Commxx const& commxx, int writer_rank)
{
    this->filename = filename;
    this->commxx = commxx;
    this->count = 0;
    this->have_file = false;
    this->iwrite_skip = write_skip;
    this->serial = serial;
    if (writer_rank == default_rank) {
        this->writer_rank = commxx.get_size() - 1;
    } else {
        this->writer_rank = writer_rank;
    }
    int idx = filename.rfind('.');
    if (idx == std::string::npos) {
        filename_base = filename;
        filename_suffix = "";
    } else {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }
    //    if (serial) {
    //        open_file();
    //    }
}

Diagnostics_write_helper::Diagnostics_write_helper(std::string const& filename,
        bool serial, int write_skip, Commxx const& commxx, int writer_rank)
{
    construct(filename, serial, write_skip, commxx, writer_rank);

}

Diagnostics_write_helper::Diagnostics_write_helper(std::string const& filename,
        bool serial, Commxx const& commxx, int writer_rank)
{
    int write_skip = 1;
    construct(filename, serial, write_skip, commxx, writer_rank);
}

Diagnostics_write_helper::Diagnostics_write_helper()
{
}

int
Diagnostics_write_helper::get_count() const
{
    return count;
}

int
Diagnostics_write_helper::get_iwrite_skip() const
{
    return iwrite_skip;
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
    return commxx.get_rank() == writer_rank;
}

int
Diagnostics_write_helper::get_writer_rank()
{
    return writer_rank;
}

H5::H5File &
Diagnostics_write_helper::get_file()
{
    if (!write_locally()) {
        throw std::runtime_error(
                "Diagnostics_write_helper::getfile() called on a non-writer rank.");
    }
    //    if (!serial) {
    open_file();
    //    }
    return *file_sptr;
}

void
Diagnostics_write_helper::finish_write()
{
    if (write_locally() && !serial) {
        file_sptr->close();
        file_sptr.reset();
        have_file = false;
    }
    ++count;
}

Diagnostics_write_helper::~Diagnostics_write_helper()
{
    if (write_locally()) {
        if (have_file) {
            file_sptr->close();
            file_sptr.reset();
            have_file = false;
        }
    }
}
