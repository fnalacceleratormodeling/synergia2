#include "synergia/bunch/diagnostics_file.h"
#include "synergia/utils/hdf5_file.h"

#include <filesystem>

namespace fs = std::filesystem;

std::string Diagnostics_file::temp_dir = "";

Diagnostics_file::Diagnostics_file(std::string const& filename,
                                   bool serial,
                                   std::shared_ptr<Commxx> const& comm)
    : file()
    , serial(serial)
    , file_count(0)
    , filename(filename)
    , filename_base()
    , filename_suffix()
    , filename_appendix()
{
    auto idx = filename.rfind('.');

    if (idx == std::string::npos) {
        filename_base = filename;
        filename_suffix = "";
    } else {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }

    file = std::make_unique<Hdf5_file>(
        get_filename(true), Hdf5_file::Flag::truncate, comm);
}

Diagnostics_file::Diagnostics_file()
    : file()
    , serial(true)
    , file_count(0)
    , filename()
    , filename_base()
    , filename_suffix()
    , filename_appendix()
{}

Diagnostics_file::~Diagnostics_file()
{
    close_file();
}

Hdf5_file&
Diagnostics_file::get_file()
{
    open_file();
    return *file;
}

std::string
Diagnostics_file::get_filename(bool use_temp_dir)
{
    std::stringstream sstream;

    if (use_temp_dir && (temp_dir != "")) {
        sstream << temp_dir;
        sstream << "/";
    }

    sstream << filename_base;

    if (!filename_appendix.empty()) { sstream << "_" << filename_appendix; }

    if (!serial) {
        sstream << "_";
        sstream << std::setw(4);
        sstream << std::setfill('0');
        sstream << file_count;
    }

    sstream << filename_suffix;

    return sstream.str();
}

void
Diagnostics_file::move_file_overwrite_if_exists(std::string const& src,
                                                std::string const& dst)
{
    fs::rename(src, dst);
}

void
Diagnostics_file::open_file()
{
    if (!file)
        file = std::make_unique<Hdf5_file>(get_filename(true),
                                           Hdf5_file::Flag::truncate);
}

void
Diagnostics_file::flush_file()
{
    if (file && (file_count % flush_period == 0)) file->flush();
}

void
Diagnostics_file::close_file()
{
    if (file) {
        file->close();

        if (!temp_dir.empty()) {
            move_file_overwrite_if_exists(get_filename(true),
                                          get_filename(false));
        }

        file.reset();
    }
}

void
Diagnostics_file::finish_write()
{
    if (serial)
        flush_file();
    else
        close_file();

    increment_count();
}
