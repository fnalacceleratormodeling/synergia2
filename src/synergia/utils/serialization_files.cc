#include "serialization_files.h"
#include "commxx.h"
#include "digits.h"

// avoid bad interaction between Boost Filesystem and clang
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>

const std::string serialization_directory("serialization");

using namespace boost::filesystem;

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
remove_directory(std::string const & name, bool parallel)
{
    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (exists(name)) {
            remove_all(name);
        }
    }
    parallel_helper.barrier();
}

void
remove_serialization_directory(bool parallel)
{
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

void
ensure_serialization_directory_exists(bool parallel)
{
    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (!is_directory(get_serialization_directory())) {
            create_directories(get_serialization_directory());
        }
    }
    parallel_helper.barrier();
}

void
rename_serialization_directory(std::string const& new_name, bool parallel)
{
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

void
symlink_serialization_directory(std::string const& existing_dir, bool parallel)
{
    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        create_symlink(existing_dir, get_serialization_directory());
    }
    parallel_helper.barrier();
}

void
unlink_serialization_directory(bool parallel)
{
    Parallel_helper parallel_helper(parallel);
    parallel_helper.barrier();
    if (parallel_helper.operate_locally()) {
        if (exists(get_serialization_directory())) {
            remove(get_serialization_directory());
        }
    }
    parallel_helper.barrier();
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
#if (defined BOOST_FILESYSTEM_VERSION) && (BOOST_FILESYSTEM_VERSION > 2)
    return get_combined_path(serialization_directory,
                             path(base_name).filename().string(), parallel);
#else
    return get_combined_path(serialization_directory,
                             path(base_name).filename(), parallel);
#endif
}

void
copy_to_serialization_directory(std::string const& file_name)
{
    copy_file_overwrite_if_exists(file_name, get_serialization_path(file_name));
}

void
copy_from_serialization_directory(std::string const& file_name)
{
    copy_file_overwrite_if_exists(get_serialization_path(file_name), file_name);
}
