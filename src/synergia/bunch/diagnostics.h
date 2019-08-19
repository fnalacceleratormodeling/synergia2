#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>
#include <map>
#include "synergia/foundation/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"

class Bunch;

/// Diagnostics is an abstract base class for bunch diagnostics classes
class Diagnostics
{

private:

    constexpr static const char * DEFAULT_WRITER_NAME = "__default_writer__";

    std::string name;
    std::string filename;
    std::string local_dir;

    bool serial;

    Bunch const * bunch;

    std::map<std::string, Diagnostics_write_helper> writers;

private:

    virtual void do_update() = 0;
    virtual void do_write() = 0;

public:

    Diagnostics( std::string const& name, bool serial,
                 std::string const& filename, 
                 std::string const& local_dir="" );

    bool have_bunch() const        { return bunch; }
    void set_bunch(Bunch const& b) { bunch = &b; }

    Bunch const& get_bunch() const
    { 
        if(!have_bunch()) throw std::runtime_error("bunch not set in diagnostics");
        return *bunch;
    }

    std::string const& get_filename()  const { return filename; }
    std::string const& get_local_dir() const { return local_dir; }

    // extra write helpers
    void delete_write_helper(std::string const & name = DEFAULT_WRITER_NAME);
    bool have_write_helper(std::string const& name = DEFAULT_WRITER_NAME) const;

    Diagnostics_write_helper &
    get_write_helper(std::string const& name = DEFAULT_WRITER_NAME);

    /// Multiple serial diagnostics can be written to a single file.
    bool is_serial() const { return serial; }

    /// Update the diagnostics
    void update() { do_update(); }

    /// Write the diagnostics to the file
    void write() { do_write(); }

    /// Update the diagnostics and write them to the file
    void update_and_write() { do_update(); do_write(); }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
};

#endif /* DIAGNOSTICS_H_ */
