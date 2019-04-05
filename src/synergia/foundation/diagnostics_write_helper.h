#ifndef DIAGNOSTICS_WRITE_HELPER_H_
#define DIAGNOSTICS_WRITE_HELPER_H_
#include <string>
#include <list>
#include "H5Cpp.h"
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
    std::string filename;
    std::string local_dir;
    bool serial;
    Commxx_sptr commxx_sptr;

    Hdf5_file_sptr file_sptr;

    bool have_file;
    int count;
    std::string filename_base, filename_suffix, filename_appendix;
    std::string
    get_filename(bool include_local_dir);
    void
    open_file();
public:
    Diagnostics_write_helper(
            std::string const & filename, 
            bool serial, 
            Commxx_sptr commxx_sptr,
            std::string const & local_dir, 
            std::string const & filename_appendix = "",
            int writer_rank = default_rank);

    // Default constructor for serialization use only
    Diagnostics_write_helper();

    /// Get the count for non-serial writers
    int
    get_count() const;

    /// Set the count for non-serial writers
    /// @param count the count
    void
    set_count(int count);

    void
    increment_count();

    bool
    write_locally();

    int
    get_writer_rank();

    Hdf5_file_sptr
    get_hdf5_file_sptr();

    void
    finish_write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    ~Diagnostics_write_helper();
};

typedef std::shared_ptr<Diagnostics_write_helper > Diagnostics_write_helper_sptr;  // syndoc:include

#endif /* DIAGNOSTICS_WRITE_HELPER_H_ */
