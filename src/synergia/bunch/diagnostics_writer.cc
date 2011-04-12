#include "diagnostics_writer.h"
#include <sstream>
#include <iomanip>

void
Diagnostics_writer::open_file_and_init()
{
    std::stringstream sstream;
    sstream << filename_base;
    if (!diagnostics_sptr->is_serial()) {
        sstream << "_";
        sstream << std::setw(4);
        sstream << std::setfill('0');
        sstream << count;
    }
    sstream << filename_suffix;
    file = H5Fcreate(sstream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
    diagnostics_sptr->init_writers(file);

}
Diagnostics_writer::Diagnostics_writer(Diagnostics *diagnostics_ptr) :
    diagnostics_ptr(diagnostics_ptr), count(0)
{
    int idx = filename.rfind('.');
    if (idx == std::string::npos) {
        filename_base = filename;
        filename_suffix = "";
    } else {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }
    if (diagnostics_ptr->is_serial()) {
        open_file_and_init();
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

void
Diagnostics_writer::write()
{
    if (!diagnostics_sptr->is_serial()) {
        open_file_and_init();
    }
    diagnostics_ptr->write();
    if (!diagnostics_sptr->is_serial()) {
        H5Fclose(file);
    }
    ++count;
}

Diagnostics_writer::~Diagnostics_writer()
{
    if (!dummy) {
        if (diagnostics_sptr->is_serial()) {
            H5Fclose(file);
        }
    }
}
