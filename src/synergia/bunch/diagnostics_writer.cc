
#include "synergia/bunch/diagnostics_writer.h"
#include "synergia/utils/hdf5_file.h"

#pragma message "TODO: replace boost::filesystem here"

Diagnostics_writer::Diagnostics_writer(
        std::string const& filename,
        std::string const& temp_dir,
        bool serial,
        Commxx const& comm )
    : file()
    , serial(serial)
    , file_count(0)
    , temp_dir(temp_dir)
    , filename(filename)
    , filename_base()
    , filename_suffix()
    , filename_appendix()
{
    auto idx = filename.rfind('.');

    if (idx == std::string::npos) 
    {
        filename_base = filename;
        filename_suffix = "";
    } 
    else 
    {
        filename_base = filename.substr(0, idx);
        filename_suffix = filename.substr(idx);
    }

    file = std::make_unique<Hdf5_file>(
            get_filename(true), Hdf5_file::truncate, comm);
}

Diagnostics_writer::Diagnostics_writer()
    : file()
    , serial(true)
    , file_count(0)
    , temp_dir()
    , filename()
    , filename_base()
    , filename_suffix()
    , filename_appendix()
{
}

Diagnostics_writer::~Diagnostics_writer()
{
    close_file();
}

Hdf5_file & Diagnostics_writer::get_file()
{
    open_file();
    return *file;
}

std::string Diagnostics_writer::get_filename(bool use_temp_dir)
{
    std::stringstream sstream;

    if (use_temp_dir && (temp_dir != "")) 
    {
        sstream << temp_dir;
        sstream << "/";
    }

    sstream << filename_base;

    if (!filename_appendix.empty()) 
    {
        sstream << "_" << filename_appendix;
    }

    if (!serial) 
    {
        sstream << "_";
        sstream << std::setw(4);
        sstream << std::setfill('0');
        sstream << file_count;
    }

    sstream << filename_suffix;

    return sstream.str();
}

void Diagnostics_writer::move_file_overwrite_if_exists(
        std::string const& src,
        std::string const& dst )
{
#if 0
    if (boost::filesystem::exists(dest)) {
        boost::filesystem::remove(dest);
    }
    boost::filesystem::copy_file(source, dest);
    boost::filesystem::remove(source);
#endif
}

void Diagnostics_writer::open_file()
{
    if (!file)
        file = std::make_unique<Hdf5_file>(
                get_filename(true), Hdf5_file::truncate);
}

void Diagnostics_writer::flush_file()
{
    if (file && (file_count % flush_period == 0))
        file->flush();
}

void Diagnostics_writer::close_file()
{
    if (file)
    {
        file->close();

        if (!temp_dir.empty())
        {
            move_file_overwrite_if_exists(
                    get_filename(true),
                    get_filename(false) );
        }

        file.reset();
    }
}

void Diagnostics_writer::finish_write()
{
    if (serial) flush_file();
    else close_file();

    increment_count();
}



