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

    if (include_local_dir && (local_dir != "")) 
    {
        sstream << local_dir;
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
        sstream << count;
    }

    sstream << filename_suffix;

    return sstream.str();
}

void
Diagnostics_write_helper::open_file()
{
    if (write_locally() && !file) 
    {
        file = std::make_unique<Hdf5_file>(
                get_filename(true).c_str(), 
                Hdf5_file::truncate );
    }
}

Diagnostics_write_helper::Diagnostics_write_helper(
        std::string const & filename,
        bool serial, 
        Commxx commxx, 
        std::string const & local_dir,
        std::string const & filename_appendix,
        int wrank ) 
    : writer_rank(wrank==default_rank ? commxx.size()-1 : wrank)
    , serial(serial)
    , count(0) 
    , commxx(commxx)
    , file()
    , local_dir(local_dir)
    , filename(filename)
    , filename_base()
    , filename_suffix()
    , filename_appendix(filename_appendix)

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
}

Hdf5_file &
Diagnostics_write_helper::get_hdf5_file()
{
    if (!write_locally()) 
    {
        throw std::runtime_error(
                "Diagnostics_write_helper::get_hdf5_file_sptr() "
                "called on a non-writer rank.");
    }

    open_file();
    return *file;
}

void
Diagnostics_write_helper::finish_write()
{
    if (write_locally() && !serial) 
    {
        file->close();

        if (local_dir != "") 
        {
            move_file_overwrite_if_exists(
                    get_filename(true), 
                    get_filename(false) );
        }

        file.reset();
    }

    ++count;

    if (write_locally() && serial) 
    {
        if (count % flush_period == 0) file->flush();
    }
}

#if 0
template<class Archive>
    void
    Diagnostics_write_helper::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & CEREAL_NVP(writer_rank)
                & CEREAL_NVP(filename)
                & CEREAL_NVP(local_dir)
                & CEREAL_NVP(serial)
                & CEREAL_NVP(commxx_sptr)
                & CEREAL_NVP(file_sptr)
                & CEREAL_NVP(have_file)
                & CEREAL_NVP(count)
                & CEREAL_NVP(filename_base)
                & CEREAL_NVP(filename_suffix)
                & CEREAL_NVP(filename_appendix);
    }

template
void
Diagnostics_write_helper::serialize<cereal::BinaryOutputArchive >(
        cereal::BinaryOutputArchive & ar, const unsigned int version);

template
void
Diagnostics_write_helper::serialize<cereal::XMLOutputArchive >(
        cereal::XMLOutputArchive & ar, const unsigned int version);

template
void
Diagnostics_write_helper::serialize<cereal::BinaryInputArchive >(
        cereal::BinaryInputArchive & ar, const unsigned int version);

template
void
Diagnostics_write_helper::serialize<cereal::XMLInputArchive >(
        cereal::XMLInputArchive & ar, const unsigned int version);
#endif

Diagnostics_write_helper::~Diagnostics_write_helper()
{
    if (write_locally() && serial && file) 
    {
        file->close();

        if (local_dir != "") 
        {
            move_file_overwrite_if_exists(
                    get_filename(true),
                    get_filename(false) );
        }
    }
}
