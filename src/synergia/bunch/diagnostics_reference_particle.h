#ifndef DIAGNOSTICS_REFERENCE_PARTICLE_H_
#define DIAGNOSTICS_REFERENCE_PARTICLE_H_

#include "synergia/bunch/diagnostics.h"

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

#endif /* DIAGNOSTICS_REFERENCE_PARTICLE_H_ */
