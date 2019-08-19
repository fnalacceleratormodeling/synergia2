#ifndef DIAGNOSTICS_WRITE_HELPER_H_
#define DIAGNOSTICS_WRITE_HELPER_H_
#include <string>
#include <list>
#include <memory>

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_file.h"

/// The Diagnostics_write_helper is a helper class for Diagnostics objects.
/// Serial Diagnostics_write_helpers write many updates to a single file.
/// Non-serial Diagnostics_write_helpers write each update to a new file.
class Diagnostics_write_helper
{
public:

    static const int default_rank = -999;
    static const int flush_period = 100;

private:

    int writer_rank;
    bool serial;
    int count;
    Commxx commxx;

    std::unique_ptr<Hdf5_file> file;

    std::string local_dir;
    std::string filename;
    std::string filename_base;
    std::string filename_suffix;
    std::string filename_appendix;

    std::string get_filename(bool include_local_dir);
    void open_file();

public:

    Diagnostics_write_helper(
            std::string const & filename, 
            bool serial, 
            Commxx commxx,
            std::string const & local_dir, 
            std::string const & filename_appendix = "",
            int writer_rank = default_rank);

    ~Diagnostics_write_helper();

    Diagnostics_write_helper(Diagnostics_write_helper &&) noexcept = default;

    /// Get the count for non-serial writers
    int get_count() const { return count; }

    /// Set the count for non-serial writers
    /// @param count the count
    void set_count(int count) { this->count = count; }
    void increment_count() { ++count; }

    bool write_locally() { return !commxx.is_null() && (writer_rank==commxx.rank()); }
    int get_writer_rank() { return writer_rank; }

    Hdf5_file& get_hdf5_file();
    void finish_write();

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

};

#endif /* DIAGNOSTICS_WRITE_HELPER_H_ */
