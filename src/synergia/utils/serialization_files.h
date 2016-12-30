#ifndef SERIALIZATION_FILES_H_
#define SERIALIZATION_FILES_H_

#include <iostream>
#include <mpi.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <cerrno>
#include <cstring>
#include <unistd.h>

#pragma clang diagnostic ignored "-Wc++11-extensions"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/archive_exception.hpp>

// copy_file_overwrite_if_exists provides portability across
// Boost Filesystem versions 2 and 3
void
copy_file_overwrite_if_exists(std::string const & source,
        std::string const & dest);

std::string
get_serialization_directory();

void
remove_directory(std::string const & name, bool parallel = true);

void
remove_serialization_directory(bool parallel = true);

void
rename_serialization_directory(std::string const& new_name,
        bool parallel = true);

void
ensure_serialization_directory_exists(bool parallel = true);

void
symlink_serialization_directory(std::string const& existing_dir,
        bool parallel = true);

void
unlink_serialization_directory(bool parallel = true);

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
    archive_save(T const& object, std::string const& filename,
            bool parallel = false)
    {
        int try_no=0;
        bool fail=true;
        while ((try_no<5) && fail){
            try {
                ensure_serialization_directory_exists(parallel);
                std::ofstream output_stream;
                output_stream.open(filename.c_str());
                int attempts=0;
                while(!output_stream.good()){
                    std::cout<<" archive_save attempt to open the file failed "<<attempts+1<<" times!"<<std::endl;
                    if (output_stream.is_open()) output_stream.close();
                    sleep(3);
                    ++attempts;
                    output_stream.open(filename.c_str(),std::ofstream::out);
                    if (attempts>5) break;
                }
                if (!output_stream.good()) {
                    std::string message("<archive>_save: unable to open ");
                    message += filename;
                    message += ": \"";
                    message += std::strerror(errno);
                    message += "\"";
                    throw std::runtime_error(message);
                }
                A output_archive(output_stream);
                int num_objects = 1;
                output_archive << BOOST_SERIALIZATION_NVP(num_objects);
                std::string object_typename(typeid(object).name());
                output_archive << BOOST_SERIALIZATION_NVP(object_typename);
                output_archive << BOOST_SERIALIZATION_NVP(object);
                output_stream.close();
                fail=false;
            }
            catch(boost::archive::archive_exception& be){
                // try again
                ++try_no;
                fail=true;
                std::cout<<" boost archive exception  has been thrown: "<<be.what()<<
                    "; trying again archive_save; trying number= "<< try_no<<std::endl<<std::flush;
                sleep(3);
            }
        }
    }

template<typename T, typename A>
    void
    archive_load(T & object, std::string const& filename)
    {
        std::ifstream input_stream(filename.c_str());
        if (!input_stream.good()) {
            std::string message("<archive>_load: unable to open ");
            message += filename;
            message += ": \"";
            message += std::strerror(errno);
            message += "\"";
            throw std::runtime_error(message);
        }
        A input_archive(input_stream);
        int num_objects;
        input_archive >> BOOST_SERIALIZATION_NVP(num_objects);
        std::string object_typename;
        input_archive >> BOOST_SERIALIZATION_NVP(object_typename);
        input_archive >> BOOST_SERIALIZATION_NVP(object);
        input_stream.close();
    }

template<typename T>
    void
    binary_save(T const& object, std::string const& filename,
            bool parallel = false)
    {
        archive_save<T, boost::archive::binary_oarchive > (object, filename,
                parallel);
    }

template<typename T>
    void
    binary_load(T & object, std::string const& filename)
    {
        archive_load<T, boost::archive::binary_iarchive > (object, filename);
    }

template<typename T>
    void
    xml_save(T const& object, std::string const& filename,
            bool parallel = false)
    {
        archive_save<T, boost::archive::xml_oarchive > (object, filename,
                parallel);
    }

template<typename T>
    void
    xml_load(T & object, std::string const& filename)
    {
        archive_load<T, boost::archive::xml_iarchive > (object, filename);
    }

#endif /* SERIALIZATION_FILES_H_ */
