#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>
#include "synergia/foundation/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"

class Bunch;

/// Diagnostics is an abstract base class for bunch diagnostics classes
class Diagnostics
{

private:

    std::string name;
    std::string filename;
    std::string local_dir;

    Bunch const * bunch;

    Diagnostics_write_helper * write_helper_ptr;
    bool have_write_helper_;

    std::map<std::string, Diagnostics_write_helper> extra_writers;

public:

    Diagnostics( std::string const& name, 
                 std::string const& filename, 
                 std::string const& local_dir="" );

    void set_bunch(Bunch const& b)
    { bunch = &b; }

    virtual std::string const&
    get_filename() const;

    virtual std::string const&
    get_local_dir() const;

    // the main write helper
    virtual void
    delete_write_helper_ptr();

    virtual Diagnostics_write_helper *
    new_write_helper_ptr();

    virtual bool
    have_write_helper() const;

    virtual Diagnostics_write_helper &
    get_write_helper();

    // extra write helpers
    void
    delete_extra_write_helper(std::string const & name);

    bool
    have_extra_write_helper(std::string const & name) const;

    Diagnostics_write_helper &
    get_extra_write_helper(std::string const & name);

    /// Multiple serial diagnostics can be written to a single file.
    virtual bool
    is_serial() const = 0;

    /// Update the diagnostics
    virtual void
    update() = 0;

    /// Write the diagnostics to the file
    virtual void
    write() = 0;

    /// Update the diagnostics and write them to the file
    virtual void
    update_and_write()
    {
        update();
        write();
    }

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Diagnostics();
};

#endif /* DIAGNOSTICS_H_ */
