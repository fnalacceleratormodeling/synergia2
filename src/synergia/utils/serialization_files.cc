#include "serialization_files.h"
#include "commxx.h"
#include <boost/filesystem.hpp>

const std::string serialization_directory("serialization");

using namespace boost::filesystem;

void
copy_file_overwrite_if_exists(std::string const & source,
        std::string const & dest)
{
    if (exists(dest)) {
        remove(dest);
    }
    copy_file(source, dest);
}

std::string
get_serialization_directory()
{
    return serialization_directory;
}

void
remove_directory(std::string const & name)
{
    Commxx commxx;
    MPI_Barrier(commxx.get());
    if (commxx.get_rank() == 0) {
        if (exists(name)) {
            remove_all(name);
        }
    }
    MPI_Barrier(commxx.get());
}

void
remove_serialization_directory()
{
    //jfa: check for link!!!!!!!
    remove_directory(get_serialization_directory());
}

void
ensure_serialization_directory_exists()
{
    Commxx commxx;
    MPI_Barrier(commxx.get());
    if (commxx.get_rank() == 0) {
        if (!is_directory(get_serialization_directory())) {
            create_directories(get_serialization_directory());
        }
    }
    MPI_Barrier(commxx.get());
}

void
rename_serialization_directory(std::string const& new_name)
{
    Commxx commxx;
    MPI_Barrier(commxx.get());
    if (commxx.get_rank() == 0) {
        if (exists(new_name)) {
            remove_all(new_name);
        }
        rename(get_serialization_directory(), new_name);
    }
    MPI_Barrier(commxx.get());
}

void
symlink_serialization_directory(std::string const& existing_dir)
{
    Commxx commxx;
    MPI_Barrier(commxx.get());
    if (commxx.get_rank() == 0) {
        create_symlink(existing_dir, get_serialization_directory());
    }
    MPI_Barrier(commxx.get());
}

void
unlink_serialization_directory()
{
    Commxx commxx;
    MPI_Barrier(commxx.get());
    if (commxx.get_rank() == 0) {
        if (exists(get_serialization_directory())) {
            remove(get_serialization_directory());
        }
    }
    MPI_Barrier(commxx.get());
}

// digits is a local function
long int
digits(long int val)
{
    long int retval = 1;
    long int base = 1;
    for (int i = 0; i < 9; ++i) {
        base *= 10;
        if (val > base) {
            retval += 1;
        }
    }
    return retval;
}

std::string
get_combined_path(std::string const& directory, std::string const& base_name,
        bool parallel)
{
    std::string full_name(base_name);
    if (parallel) {
        Commxx commxx;
        std::stringstream sstream;
        sstream << std::setw(digits(commxx.get_size()));
        sstream << std::setfill('0');
        sstream << commxx.get_rank();
        sstream << "_";
        sstream << base_name;
        full_name = sstream.str();
    }
    return directory + "/" + full_name;
}

std::string
get_serialization_path(std::string const& base_name, bool parallel)
{
    return get_combined_path(serialization_directory, base_name, parallel);
}

void
copy_to_serialization_directory(std::string const& file_name)
{
    //    ensure_serialization_directory_exists();
    copy_file_overwrite_if_exists(file_name, get_serialization_path(file_name));
}

void
copy_from_serialization_directory(std::string const& file_name)
{
//    ensure_serialization_directory_exists();
    copy_file_overwrite_if_exists(get_serialization_path(file_name), file_name);
}
