#ifndef DIAGNOSTICS_BULK_TRACK_H_
#define DIAGNOSTICS_BULK_TRACK_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_bulk_track records the phase space coordinates of a multiple
/// particles.
/// Particles will only be tracked if they stay on the same processor.
/// Lost particles that are somehow restored or particles not available when
/// the first update is called will also not be tracked.
class Diagnostics_bulk_track : public Diagnostics
{

private:

    int total_num_tracks, local_num_tracks;
    int offset, local_offset;
    bool setup;

    // used between update and write
    double s_n;
    int repetition;
    double s;
    double pz;

    double ref_charge;
    double ref_mass;
    double ref_pz;

    karray2d_row track_coords;

private:

    void do_update(Bunch const& bunch) override;
    void do_reduce(Commxx comm, int writer_rank) override { }
    void do_first_write(Hdf5_file & file) override;
    void do_write(Hdf5_file & file) override;

    friend class cereal::access;

    template<class AR>
    void serialize(AR & ar)
    { 
        ar(cereal::base_class<Diagnostics>(this));
        ar(total_num_tracks);
        ar(local_num_tracks);
        ar(offset);
        ar(local_offset);
        ar(setup);
    }

public:

    /// Create an empty Diagnostics_bulk_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param num_tracks the number of local particles to track
    /// @param offset id offset for first particle to track
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_bulk_track(int num_tracks = 0, int offset = 0);

};

CEREAL_REGISTER_TYPE(Diagnostics_bulk_track)

#endif /* DIAGNOSTICS_BULK_TRACK_H_ */
