#ifndef SERIALIZATION_FILES_H_
#define SERIALIZATION_FILES_H_

#include <string>

std::string
get_serialization_directory();

void
ensure_serialization_directory_exists();

std::string
get_serialization_path(std::string const& base_name);

void
copy_to_serialization_directory(std::string const& file_name);

void
copy_from_serialization_directory(std::string const& file_name);

#endif /* SERIALIZATION_FILES_H_ */
