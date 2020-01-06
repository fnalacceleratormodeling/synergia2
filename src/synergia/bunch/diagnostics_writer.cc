
#include "synergia/bunch/diagnostics_writer.h"
#include "synergia/utils/hdf5_file.h"

Diagnostics_writer::Diagnostics_writer(
        std::string const& filename,
        std::string const& local_dir)
{

}

Diagnostics_writer::Diagnostics_writer()
{
}

Diagnostics_writer::~Diagnostics_writer()
{
    close_file();
}

Hdf5_file & Diagnostics_writer::get_file()
{
    return *file;
}

std::string Diagnostics_writer::get_filename(bool include_local_dir)
{
    return "";
}

void Diagnostics_writer::move_file_overwrite_if_exists()
{
}

bool Diagnostics_writer::write_locally()
{
    return true;
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

        if (!local_dir.empty())
        {
            move_file_overwrite_if_exists();
        }

        file.reset();
    }
}

void Diagnostics_writer::finish_write(bool serial)
{
    if (serial) flush_file();
    else close_file();
}



