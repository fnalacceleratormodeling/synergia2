#ifndef SERIALIZATION_FILES_H_
#define SERIALIZATION_FILES_H_

#include <fstream>
#include <stdexcept>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

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

template<typename T, typename A>
    void
    archive_save(T const& object, std::string const& filename)
    {
        std::ofstream output_stream(filename.c_str());
        if (!output_stream.good()) {
            std::string message("<archive>_save: unable to open ");
            message += filename;
            throw std::runtime_error(message);
        }
        ensure_serialization_directory_exists();
        A output_archive(output_stream);
        output_archive << BOOST_SERIALIZATION_NVP(object);
        output_stream.close();
    }

template<typename T, typename A>
    void
    archive_load(T & object, std::string const& filename)
    {
        std::ifstream input_stream(filename.c_str());
        if (!input_stream.good()) {
            std::string message("<archive>_load: unable to open ");
            message += filename;
            throw std::runtime_error(message);
        }
        A input_archive(input_stream);
        input_archive >> BOOST_SERIALIZATION_NVP(object);
        input_stream.close();
    }

template<typename T>
    void
    binary_save(T const& object, std::string const& filename)
    {
        archive_save<T, boost::archive::binary_oarchive > (object, filename);
    }

template<typename T>
    void
    binary_load(T & object, std::string const& filename)
    {
        archive_load<T, boost::archive::binary_iarchive > (object, filename);
    }

template<typename T>
    void
    text_save(T const& object, std::string const& filename)
    {
        archive_save<T, boost::archive::text_oarchive > (object, filename);
    }

template<typename T>
    void
    text_load(T & object, std::string const& filename)
    {
        archive_load<T, boost::archive::text_iarchive > (object, filename);
    }

template<typename T>
    void
    xml_save(T const& object, std::string const& filename)
    {
        archive_save<T, boost::archive::xml_oarchive > (object, filename);
    }

template<typename T>
    void
    xml_load(T & object, std::string const& filename)
    {
        archive_load<T, boost::archive::xml_iarchive > (object, filename);
    }

#endif /* SERIALIZATION_FILES_H_ */
