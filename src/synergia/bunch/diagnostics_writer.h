#ifndef DIAGNOSTICS_WRITER_H_
#define DIAGNOSTICS_WRITER_H_
#include "synergia/bunch/diagnostics.h"
#include <string>
#include <list>

/// The Diagnostics_writer writes a Diagnostics object to a file.
/// Serial Diagnostics_writers write many updates to a single file.
/// Non-serial Diagnostics_writers write each update to a new file.
//class Diagnostics_writer
//{
//private:
//    Diagnostics_sptr diagnostics_sptr;
//    hid_t file;
//    bool dummy;
//    int count;
//    std::string filename_base, filename_suffix;
//    void
//    open_file_and_init();
//public:
//    /// Construct Diagnostics_writer
//    /// @param filename the filename to which to write
//    /// @param diagnostics_sptr the diagnostics object to write
//    Diagnostics_writer(std::string const& filename,
//            Diagnostics_sptr diagnostics_sptr);
//
//    /// Construct a dummy Diagnostics_writer
//    Diagnostics_writer();
//
//    /// Whether the Diagnostics_writer is really a dummy
//    bool
//    is_dummy() const;
//
//    /// Get the Diagnostics_sptr to be used
//    Diagnostics_sptr
//    get_diagnostics_sptr();
//
//    /// Get the count for non-serial writers
//    int
//    get_count() const;
//
//    /// Set the count for non-serial writers
//    /// @param count the count
//    void
//    set_count(int count);
//
//    /// Write the Diagnostics to the file
//    void
//    write();
//
//    /// Update the Diagnostics with the state of a Bunch and write to the file
//    /// @param bunch the Bunch to be passed to the Diagnostics
//    void
//    update_and_write(Bunch const& bunch);
//
//    ~Diagnostics_writer();
//};
//
//typedef boost::shared_ptr<Diagnostics_writer > Diagnostics_writer_sptr;
//
///// A pre-defined dummy Diagnostics_writer
//Diagnostics_writer
//no_diagnostics();

/// Multi_diagnostics contains a list of Diagnostics_sptrs
class Multi_diagnostics
{
private:
    std::list<Diagnostics_sptr > diagnostics;
public:
    /// Construct an empty list of Diagnostics_writer_sptrs
    Multi_diagnostics();

    /// Append a Diagnostics_sptr to the list
    /// @param diagnostics_sptr the Diagnostics_sptr
    void append(Diagnostics_sptr diagnostics_sptr);

    /// The same as append -- included for notational consistency with C++
    /// @param diagnostics_sptr the Diagnostics_sptr
    void push_back(Diagnostics_sptr diagnostics_sptr);

    /// A convenience definition of the list iterator. Not relevant for
    /// Python.
    typedef std::list<Diagnostics_sptr >::iterator iterator;

    /// A convenience definition of the list iterator begin. Not relevant for
    /// Python.
    iterator begin();

    /// A convenience definition of the list iterator end. Not relevant for
    /// Python.
    iterator end();
};

///// A pre-defined dummy Multi_diagnostics
Multi_diagnostics
no_diagnostics();


#endif /* DIAGNOSTICS_WRITER_H_ */
