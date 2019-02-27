#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"

/// Diagnostics is an abstract base class for bunch diagnostics classes
class Diagnostics
{
private:
    std::string name;
    std::string filename;
    std::string local_dir;
    Bunch_sptr bunch_sptr;
    bool have_bunch_;
    Diagnostics_write_helper * write_helper_ptr;
    bool have_write_helper_;

public:

    Diagnostics(std::string const& name, std::string const& filename, std::string const& local_dir="");

    Diagnostics(Diagnostics const&) = delete;
    Diagnostics& operator=(Diagnostics const&) = delete;

    // Default constructor for serialization use only
    Diagnostics();

    virtual std::string const&
    get_filename() const;

    virtual std::string const&
    get_local_dir() const;

    virtual void
    set_bunch_sptr(Bunch_sptr bunch_sptr);

    virtual bool
    have_bunch() const;

    virtual void
    delete_write_helper_ptr();

    virtual Diagnostics_write_helper *
    new_write_helper_ptr();

    virtual bool
    have_write_helper() const;

    virtual Diagnostics_write_helper &
    get_write_helper();

    Bunch &
    get_bunch();
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
    update_and_write();

    template<class Archive>
     void
    serialize(Archive & ar, const unsigned int version);

    virtual
    ~Diagnostics();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics)

typedef boost::shared_ptr<Diagnostics > Diagnostics_sptr; // syndoc:include
typedef std::list<Diagnostics_sptr > Diagnosticss;
typedef std::vector<Diagnosticss > Train_diagnosticss;

#endif /* DIAGNOSTICS_H_ */
