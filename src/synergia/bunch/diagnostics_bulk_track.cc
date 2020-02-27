#include <string>
#include <iostream>
#include "diagnostics_bulk_track.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/simple_timer.h"


Diagnostics_bulk_track::Diagnostics_bulk_track(
        int num_tracks, 
        int offset )
    : Diagnostics("diagnostis_bulk_track", true)
    , total_num_tracks(num_tracks)
    , local_num_tracks(0)
    , offset(offset)
    , local_offset(0)
    , setup(false)
    , track_coords("local_coords", 0, 0)
{
}

void
Diagnostics_bulk_track::do_update(Bunch const& bunch)
{
    scoped_simple_timer("diag_bulk_track_update");

    auto const& ref = bunch.get_reference_particle();

    if (!setup)
    {
        ref_charge = ref.get_charge();
        ref_mass   = ref.get_four_momentum().get_mass();
        ref_pz     = ref.get_four_momentum().get_momentum();

        auto const& comm = bunch.get_comm();

        local_num_tracks = decompose_1d_local(comm, total_num_tracks);
        local_offset = decompose_1d_local(comm, offset);

        if (local_num_tracks + local_offset > bunch.get_local_num_slots()) 
            local_num_tracks = bunch.get_local_num_slots() - local_offset;

        setup = true;
    }

    pz         = ref.get_momentum();
    s          = ref.get_s();
    s_n        = ref.get_s_n();
    repetition = ref.get_repetition();

    track_coords = bunch.get_particles_in_range(local_offset, local_num_tracks);
}

void
Diagnostics_bulk_track::do_first_write(Hdf5_file& file)
{
    file.write("charge", ref_charge);
    file.write("mass", ref_mass);
    file.write("pz", ref_pz);
}

void
Diagnostics_bulk_track::do_write(Hdf5_file& file)
{   
    scoped_simple_timer("diag_bulk_track_write");

    // write serial
    file.append_single("track_pz", pz);
    file.append_single("track_s", s);
    file.append_single("track_s_n", s_n);
    file.append_single("track_repetition", repetition);

    // write collective from all ranks
    file.append_collective("track_coords", track_coords);
}


