#ifndef DIAGNOSTICS_BULK_TRACK_H_
#define DIAGNOSTICS_BULK_TRACK_H_

#include "synergia/bunch/diagnostics.h"
#include "synergia/utils/logger.h"

/// Diagnostics_bulk_track records the phase space coordinates of a multiple
/// particles.
/// Particles will only be tracked if they stay on the same processor.
/// Lost particles that are somehow restored or particles not available when
/// the first update is called will also not be tracked.
class Diagnostics_bulk_track : public Diagnostics
{
public:

    constexpr static const char* diag_type = "diagnostics_bulk_track";
    constexpr static const bool  diag_write_serial = true;

private:

    struct Track_status 
    {
      int index;
      int particle_id;

      template<class Archive>
      void serialize(Archive & ar, const unsigned int version);
    };

    int total_num_tracks, local_num_tracks;
    int offset, local_offset;
    bool first_search;
    bool first_write;
    std::vector<Track_status> diag_track_status;

    std::vector<int> num_tracks;
    std::vector<int> track_displs;

    double s_n;
    int repetition;
    double s;
    double pz;

    karray2d_row track_coords;
    karray2d_row all_coords;

private:

    void receive_other_local_coords(std::vector<int > const& local_nums);
    void send_local_coords();

    void do_update(Bunch const& bunch) override;
    void do_write (Bunch const& bunch) override;

public:

    /// Create an empty Diagnostics_bulk_track object
    /// @param bunch_sptr the Bunch
    /// @param filename the base name for file to write to (base names will have
    ///        a numerical index inserted
    /// @param num_tracks the number of local particles to track
    /// @param offset id offset for first particle to track
    /// @param local_dir local directory to use for temporary scratch
    Diagnostics_bulk_track(
            int num_tracks,
            int offset,
            std::string const& filename,
			std::string const& local_dir = "" );

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

};

#endif /* DIAGNOSTICS_BULK_TRACK_H_ */
