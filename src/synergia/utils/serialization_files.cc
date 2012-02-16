#include "serialization_files.h"
#include <boost/filesystem.hpp>

const std::string serialization_directory("serialization");

std::string
get_serialization_directory()
{
    return serialization_directory;
}

void
ensure_serialization_directory_exists()
{
    if (!boost::filesystem::is_directory(get_serialization_directory())) {
        boost::filesystem::create_directories(get_serialization_directory());
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
    boost::filesystem::copy_file(file_name, get_serialization_path(file_name),
            boost::filesystem::copy_option::overwrite_if_exists);
}

void
copy_from_serialization_directory(std::string const& file_name)
{
    ensure_serialization_directory_exists();
    boost::filesystem::copy_file(get_serialization_path(file_name), file_name,
            boost::filesystem::copy_option::overwrite_if_exists);
}
