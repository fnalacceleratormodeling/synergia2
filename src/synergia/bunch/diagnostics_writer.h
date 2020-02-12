#ifndef DIAGNOSTICS_WRITER_H
#define DIAGNOSTICS_WRITER_H

#include <memory>
#include <cereal/types/memory.hpp>


class Commxx;
class Hdf5_file;

class Diagnostics_writer
{
private:

    std::unique_ptr<Hdf5_file> file;

    bool serial;
    int file_count;

    std::string temp_dir;
    std::string filename;
    std::string filename_base;
    std::string filename_suffix;
    std::string filename_appendix;

    static const int flush_period = 100;

public:

    Diagnostics_writer( 
            std::string const& filename,
            std::string const& temp_dir,
            bool serial,
            Commxx const& comm);

    ~Diagnostics_writer();

    Diagnostics_writer(Diagnostics_writer const&) = delete;
    Diagnostics_writer(Diagnostics_writer &&) noexcept = default;

    Hdf5_file & get_file();

    // count for non-serial writer
    int  get_count() const     { return file_count; }
    void set_count(int count)  { file_count = count; }
    void increment_count()     { ++file_count; }

    void finish_write();

    void open_file();
    void flush_file();
    void close_file();

private:

    std::string get_filename(
            bool use_temp_dir);

    void move_file_overwrite_if_exists(
            std::string const& src, 
            std::string const& dst );

    friend class Diagnostics_worker;
    friend class cereal::access;

    Diagnostics_writer();

    template<class AR>
    void serialize(AR & ar)
    {
        ar(CEREAL_NVP(file));
        ar(CEREAL_NVP(serial));
        ar(CEREAL_NVP(file_count));
        ar(CEREAL_NVP(temp_dir));
        ar(CEREAL_NVP(filename));
        ar(CEREAL_NVP(filename_base));
        ar(CEREAL_NVP(filename_suffix));
        ar(CEREAL_NVP(filename_appendix));
    }
};

#endif
