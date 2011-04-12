#ifndef DIAGNOSTICS_WRITER_H_
#define DIAGNOSTICS_WRITER_H_
#include <string>
#include <list>
#include "hdf5.h"
#include <boost/shared_ptr.hpp>

class Diagnostics;

/// The Diagnostics_writer is a helper class for Diagnostics objects.
/// Serial Diagnostics_writers write many updates to a single file.
/// Non-serial Diagnostics_writers write each update to a new file.
class Diagnostics_writer
{
private:
    Diagnostics *diagnostics_ptr;
    hid_t file;
    int count;
    std::string filename_base, filename_suffix;
    void
    open_file_and_init();
public:
    /// Construct Diagnostics_writer
    Diagnostics_writer(Diagnostics *diagnostics_ptr);

    /// Get the count for non-serial writers
    int
    get_count() const;

    /// Set the count for non-serial writers
    /// @param count the count
    void
    set_count(int count);

    /// Write the Diagnostics to the file
    void
    write();

    ~Diagnostics_writer();
};

typedef boost::shared_ptr<Diagnostics_writer > Diagnostics_writer_sptr;

#endif /* DIAGNOSTICS_WRITER_H_ */
