#ifndef DIAGNOSTICS_BULK_SPECTATOR_TRACK_H_
#define DIAGNOSTICS_BULK_SPECTATOR_TRACK_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/logger.h"

/// Diagnostics_bulk_spectator_track records the phase space coordinates of a multiple
/// particles.
/// Particles will only be tracked if they stay on the same processor.
/// Lost particles that are somehow restored or particles not available when
/// the first update is called will also not be tracked.
class Diagnostics_bulk_spectator_track : public Diagnostics
{
public:
    static const char name[];
private:
    struct Track_status {
      bool found;
      int last_index;
      int particle_id;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    };

    int num_tracks, local_num_tracks;
    int offset, local_offset;
    bool have_writers;
    bool first_search;
    std::vector<Track_status > diag_track_status;
    double s_n;
    Hdf5_serial_writer<double > * writer_s_n;
    int repetition;
    Hdf5_serial_writer<int > * writer_repetition;
    double s;
    Hdf5_serial_writer<double > * writer_s;
    double pz;
    Hdf5_serial_writer<double > * writer_pz;
    MArray2d track_coords;
    Hdf5_serial_writer<MArray2d_ref > * writer_coords;
    virtual void
    init_writers(Hdf5_file_sptr file_sptr);
    void
    receive_other_local_coords(std::vector<int > const& local_nums);
    void
    send_local_coords();
public:
    /// Create an empty Diagnostics_bulk_spectator_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param num_tracks the number of local particles to track
    /// @param offset id offset for first particle to track
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_bulk_spectator_track(std::string const& filename,
			   int num_tracks, int offset=0, std::string const& local_dir = "");

    // Default constructor for serialization use only
    Diagnostics_bulk_spectator_track();

//    virtual Diagnostics_write_helper *
//    new_write_helper_ptr();

    /// Multiple serial diagnostics can be written to a single file.
    /// The Diagnostics_bulk_spectator_track class is serial.
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
    ~Diagnostics_bulk_spectator_track();
};
BOOST_CLASS_EXPORT_KEY(Diagnostics_bulk_spectator_track)
typedef boost::shared_ptr<Diagnostics_bulk_spectator_track > Diagnostics_bulk_spectator_track_sptr; // syndoc:include

#endif /* DIAGNOSTICS_BULK_SPECTATOR_TRACK_H_ */
