#ifndef SERIALIZATION_FILES_H_
#define SERIALIZATION_FILES_H_

#include <string>

// copy_file_overwrite_if_exists provides portability across
// Boost Filesystem versions 2 and 3
void
copy_file_overwrite_if_exists(std::string const & source,
        std::string const & dest);

std::string
get_serialization_directory();

void
remove_directory(std::string const & name);

void
remove_serialization_directory();

void
rename_serialization_directory(std::string const& new_name);

void
ensure_serialization_directory_exists();

void
symlink_serialization_directory(std::string const& existing_dir);

void
unlink_serialization_directory();

std::string
get_combined_path(std::string const& directory, std::string const& base_name,
        bool parallel = true);

std::string
get_serialization_path(std::string const& base_name, bool parallel = true);

void
copy_to_serialization_directory(std::string const& file_name);

void
copy_from_serialization_directory(std::string const& file_name);

#endif /* SERIALIZATION_FILES_H_ */
