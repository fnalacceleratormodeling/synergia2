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
    Bunch_sptr bunch_sptr;
    bool have_bunch_;
    Diagnostics_write_helper * write_helper_ptr;
    bool have_write_helper_;

public:

    Diagnostics(std::string const& name, std::string const& filename);

    // Default constructor for serialization use only
    Diagnostics();

    virtual std::string const&
    get_filename() const;

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
    update_and_write()
    {
        update();
        write();
    }
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(name);
            ar & BOOST_SERIALIZATION_NVP(filename);
            ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
            ar & BOOST_SERIALIZATION_NVP(have_bunch_);
            ar & BOOST_SERIALIZATION_NVP(write_helper_ptr);
            ar & BOOST_SERIALIZATION_NVP(have_write_helper_);
        }
    virtual
    ~Diagnostics();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics)
typedef boost::shared_ptr<Diagnostics > Diagnostics_sptr;



class Diagnostics_reference_particle : public Diagnostics
{
public:
    static const char name[];
private:
    bool have_writers;
    Hdf5_serial_writer<double > * writer_beta;
    Hdf5_serial_writer<double > * writer_gamma;
    Hdf5_serial_writer<MArray1d_ref > * writer_state;
    Hdf5_serial_writer<double > * writer_s;

public:
    /// Create a Diagnostics_reference_particle object
    /// @param bunch the Bunch
    /// @param filename filename for output
    Diagnostics_reference_particle(
            std::string const& filename);

    // Default constructor for serialization use only
    Diagnostics_reference_particle();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_reference_particle class is serial.
    virtual bool
    is_serial() const;

    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                    & BOOST_SERIALIZATION_NVP(have_writers)
                    & BOOST_SERIALIZATION_NVP(writer_beta)
                    & BOOST_SERIALIZATION_NVP(writer_gamma)
                    & BOOST_SERIALIZATION_NVP(writer_state)
                    & BOOST_SERIALIZATION_NVP(writer_s);
        }

    virtual
    ~Diagnostics_reference_particle();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_reference_particle)
typedef boost::shared_ptr<Diagnostics_reference_particle >
        Diagnostics_reference_particle_sptr;

#endif /* DIAGNOSTICS_H_ */
