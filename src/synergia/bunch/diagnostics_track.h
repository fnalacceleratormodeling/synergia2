#ifndef DIAGNOSTICS_TRACK_H_
#define DIAGNOSTICS_TRACK_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_track records the phase space coordinates of a single particle.
/// Particles will only be tracked if they stay on the same processor.
/// Lost particles that are somehow restored or particles not available when
/// the first update is called will also not be tracked.
class Diagnostics_track : public Diagnostics
{
public:
    static const char name[];
private:
    bool have_writers;
    bool found;
    bool first_search;
    int last_index;
    int particle_id;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double trajectory_length;
    Hdf5_serial_writer<double > * writer_trajectory_length;
    MArray1d coords;
    Hdf5_serial_writer<MArray1d_ref > * writer_coords;
    virtual void
    init_writers(Hdf5_file_sptr file_sptr);

public:
    /// Create an empty Diagnostics_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param particle_id the particle ID to track
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_track(std::string const& filename,
            int particle_id, std::string const& local_dir="");

    // Default constructor for serialization use only
    Diagnostics_track();

    virtual Diagnostics_write_helper *
    new_write_helper_ptr();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_track class is serial.
    virtual bool
    is_serial() const;

    /// Update the diagnostics
    virtual void
    update();

    virtual void
    write();

    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);

    virtual
    ~Diagnostics_track();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_track)
typedef boost::shared_ptr<Diagnostics_track > Diagnostics_track_sptr; // syndoc:include

#endif /* DIAGNOSTICS_TRACK_H_ */
