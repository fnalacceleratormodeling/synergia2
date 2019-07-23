#ifndef UTILS_CEREAL_FILES_H_
#define UTILS_CEREAL_FILES_H_

#include <iostream>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <cerrno>
#include <cstring>
#include <unistd.h>
#include <iomanip>

#pragma clang diagnostic ignored "-Wc++11-extensions"

#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

#include "commxx.h"
#include "digits.h"

// avoid bad interaction between Boost Filesystem and clang
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>

constexpr const char * serialization_directory = "serialization";

class Parallel_helper
{
private:
    bool parallel;
public:
    Parallel_helper(bool parallel) :
        parallel(parallel)
    {
    }
    void
    barrier()
    {
        if (parallel) {
            MPI_Barrier(Commxx());
        }
    }
    bool
    operate_locally()
    {
        if (parallel) {
            return Commxx().get_rank() == 0;
        } else {
            return true;
        }
    }
};



// copy_file_overwrite_if_exists provides portability across
// Boost Filesystem versions 2 and 3
inline void
copy_file_overwrite_if_exists(std::string const & source, std::string const & dest)
{
    using namespace boost::filesystem;

    if (exists(dest)) {
        remove(dest);
    }
    copy_file(source, dest);
}

inline std::string
get_serialization_directory()
{
    return serialization_directory;
}

inline void
remove_directory(std::string const & name, bool parallel = true)
{
    using namespace boost::filesystem;

    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (exists(name)) {
            remove_all(name);
        }
    }
    parallel_helper.barrier();
}

inline void
remove_serialization_directory(bool parallel = true)
{
    using namespace boost::filesystem;

    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (is_symlink(get_serialization_directory())) {
            remove(get_serialization_directory());
        } else {
            if (exists(get_serialization_directory())) {
                remove_all(get_serialization_directory());
            }
        }
    }
    parallel_helper.barrier();
}

inline void
rename_serialization_directory(std::string const& new_name, bool parallel = true)
{
    using namespace boost::filesystem;

    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (exists(new_name)) {
            remove_all(new_name);
        }
        rename(get_serialization_directory(), new_name);
    }
    parallel_helper.barrier();
}

inline void
ensure_serialization_directory_exists(bool parallel = true)
{
    using namespace boost::filesystem;

    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (!is_directory(get_serialization_directory())) {
            create_directories(get_serialization_directory());
        }
    }
    parallel_helper.barrier();
}



inline void
symlink_serialization_directory(std::string const& existing_dir, bool parallel = true)
{
    using namespace boost::filesystem;

    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        create_symlink(existing_dir, get_serialization_directory());
    }
    parallel_helper.barrier();
}

inline void
unlink_serialization_directory(bool parallel = true)
{
    using namespace boost::filesystem;

    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (exists(get_serialization_directory())) {
            remove(get_serialization_directory());
        }
    }
    parallel_helper.barrier();
}


inline std::string
get_combined_path(std::string const& directory, 
        std::string const& base_name, 
        bool parallel = true)
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


inline std::string
get_serialization_path(std::string const& base_name, bool parallel = true)
{
    using namespace boost::filesystem;

#if (defined BOOST_FILESYSTEM_VERSION) && (BOOST_FILESYSTEM_VERSION > 2)
    return get_combined_path(serialization_directory,
                             path(base_name).filename().string(), parallel);
#else
    return get_combined_path(serialization_directory,
                             path(base_name).filename(), parallel);
#endif
}

inline void
copy_to_serialization_directory(std::string const& file_name)
{
    copy_file_overwrite_if_exists(file_name, get_serialization_path(file_name));
}


inline void
copy_from_serialization_directory(std::string const& file_name)
{
    copy_file_overwrite_if_exists(get_serialization_path(file_name), file_name);
}

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
                output_archive << CEREAL_NVP(num_objects);
                std::string object_typename(typeid(object).name());
                output_archive << CEREAL_NVP(object_typename);
                output_archive << CEREAL_NVP(object);
                output_stream.close();
                fail=false;
            }
            //catch(boost::archive::archive_exception& be){
            catch(...){
                // try again
                ++try_no;
                fail=true;
                std::cout<<" boost archive exception  has been thrown: "<<//be.what()<<
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
        input_archive >> CEREAL_NVP(num_objects);
        std::string object_typename;
        input_archive >> CEREAL_NVP(object_typename);
        input_archive >> CEREAL_NVP(object);
        input_stream.close();
    }

template<typename T>
    void
    binary_save(T const& object, std::string const& filename, bool parallel = false)
    {
        archive_save<T, cereal::BinaryOutputArchive> (object, filename, parallel);
    }

template<typename T>
    void
    binary_load(T & object, std::string const& filename)
    {
        archive_load<T, cereal::BinaryInputArchive> (object, filename);
    }

template<typename T>
    void
    xml_save(T const& object, std::string const& filename, bool parallel = false)
    {
        archive_save<T, cereal::XMLOutputArchive> (object, filename, parallel);
    }

template<typename T>
    void
    xml_load(T & object, std::string const& filename)
    {
        archive_load<T, cereal::XMLInputArchive> (object, filename);
    }

#endif /* SERIALIZATION_FILES_H_ */
