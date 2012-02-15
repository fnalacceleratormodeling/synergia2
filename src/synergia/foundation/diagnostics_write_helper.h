#ifndef DIAGNOSTICS_WRITE_HELPER_H_
#define DIAGNOSTICS_WRITE_HELPER_H_
#include <string>
#include <list>
#include "H5Cpp.h"
#include <boost/shared_ptr.hpp>

#include "synergia/utils/commxx.h"
#include "synergia/utils/hdf5_file.h"

/// The Diagnostics_write_helper is a helper class for Diagnostics objects.
/// Serial Diagnostics_write_helpers write many updates to a single file.
/// Non-serial Diagnostics_write_helpers write each update to a new file.
class Diagnostics_write_helper
{
public:
    static const int default_rank = -999;
private:
    int writer_rank;
    std::string filename;
    bool serial;
    Commxx commxx;
    Hdf5_file_sptr file_sptr;
    bool have_file;
    int count;
    int iwrite_skip; // skip writing iwrite_skip turns or steps
    std::string filename_base, filename_suffix;
    std::string
    get_filename();
    void
    open_file();
    void
    construct(std::string const& filename,
        bool serial, int write_skip, Commxx const& commxx, int writer_rank);
public:
    /// Construct Diagnostics_write_helper
    Diagnostics_write_helper(std::string const& filename, bool serial, int write_skip, Commxx const& commxx,
            int writer_rank = default_rank);

    Diagnostics_write_helper(std::string const& filename, bool serial, Commxx const& commxx,
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

    void increment_count();

    int get_iwrite_skip()const;

    bool
    write_locally();

    int
    get_writer_rank();

    Hdf5_file_sptr
    get_hdf5_file_sptr();

    void
    finish_write();

    ~Diagnostics_write_helper();
};

typedef boost::shared_ptr<Diagnostics_write_helper > Diagnostics_write_helper_sptr;

#endif /* DIAGNOSTICS_WRITE_HELPER_H_ */
