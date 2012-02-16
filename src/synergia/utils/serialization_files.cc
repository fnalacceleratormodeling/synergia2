#include "serialization_files.h"
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
ensure_serialization_directory_exists()
{
    if (!is_directory(get_serialization_directory())) {
        create_directories(get_serialization_directory());
    }
}

std::string
get_serialization_path(std::string const& base_name)
{
    return serialization_directory + "/" + base_name;
}

void
copy_to_serialization_directory(std::string const& file_name)
{
    ensure_serialization_directory_exists();
    copy_file_overwrite_if_exists(file_name, get_serialization_path(file_name));
}

void
copy_from_serialization_directory(std::string const& file_name)
{
    ensure_serialization_directory_exists();
    copy_file_overwrite_if_exists(get_serialization_path(file_name), file_name);
}
