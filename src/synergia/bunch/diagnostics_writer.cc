#include "diagnostics_writer.h"
#include <sstream>
#include <iomanip>

void
Diagnostics_writer::open_file()
{
    if (!have_file) {
        std::stringstream sstream;
        sstream << filename_base;
        if (!serial) {
            sstream << "_";
            sstream << std::setw(4);
            sstream << std::setfill('0');
            sstream << count;
        }
        sstream << filename_suffix;
        file = H5Fcreate(sstream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                H5P_DEFAULT);
        have_file = true;
    }
}

Diagnostics_writer::Diagnostics_writer(std::string const& filename,
        bool serial, Commxx const& commxx) :
    filename(filename), serial(serial), commxx(commxx), count(0), have_file(
            false)
{
    int idx = filename.rfind('.');
    if (idx == std::string::npos) {
        filename_base = filename;
        filename_suffix = "";
    } else {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }
    if (serial) {
        open_file();
    }
}

int
Diagnostics_writer::get_count() const
{
    return count;
}

void
Diagnostics_writer::set_count(int count)
{
    this->count = count;
}

bool
Diagnostics_writer::write_locally()
{
    return commxx.get_rank() == 0;
}

hid_t &
Diagnostics_writer::get_hdf5_file()
{
    if (!serial) {
        open_file();
    }
    return file;
}

void
Diagnostics_writer::finish_write()
{
    if (!serial) {
        H5Fclose(file);
        have_file = false;
    }
    ++count;
}

Diagnostics_writer::~Diagnostics_writer()
{
    if (serial) {
        H5Fclose(file);
        have_file = false;
    }
}
