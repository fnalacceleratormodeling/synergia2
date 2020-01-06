#ifndef DIAGNOSTICS_WRITER_H
#define DIAGNOSTICS_WRITER_H

#include <memory>

#include <cereal/types/memory.hpp>

class Hdf5_file;

class Diagnostics_writer
{
private:

    std::unique_ptr<Hdf5_file> file;

    int file_count;
    int wrank;

    std::string local_dir;
    std::string filename;
    std::string filename_base;
    std::string filename_suffix;
    std::string filename_appendix;

    static const int flush_period = 100;

public:

    Diagnostics_writer( 
            std::string const& filename,
            std::string const& local_dir = "" );

    Diagnostics_writer();
    ~Diagnostics_writer();

    Diagnostics_writer(Diagnostics_writer const&) = delete;
    Diagnostics_writer(Diagnostics_writer &&) noexcept = default;

    Hdf5_file & get_file();
    int writer_rank() const { return wrank; }

    // count for non-serial writer
    int  get_count() const { return file_count; }
    void set_count(int count) { file_count = count; }
    void increment_count() { ++file_count; }

    bool write_locally();
    void finish_write(bool serial);

    void open_file();
    void flush_file();
    void close_file();

private:

    std::string get_filename(bool include_local_dir);
    void move_file_overwrite_if_exists();

    friend class cereal::access;

    template<class AR>
    void serialize(AR & ar)
    {
    }
};

#endif
