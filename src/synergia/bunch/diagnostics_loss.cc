#include "synergia/bunch/diagnostics_loss.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/simple_timer.h"


void Diagnostics_loss::do_update(Bunch const& bunch)
{ 
    scoped_simple_timer timer("diag_loss_update");

    auto ref = bunch.get_reference_particle();
    repetition      = ref.get_repetition();
    s_ref_particle  = ref.get_s();
    sn_ref_particle = ref.get_s_n();
    bucket_index    = bunch.get_bucket_index();

    coords = bunch.get_particles_last_discarded();
}

void Diagnostics_loss::do_first_write(Hdf5_file & file)
{ 
    file.append("bucket_index", bucket_index);
}

void Diagnostics_loss::do_write(Hdf5_file & file)
{ 
    scoped_simple_timer timer("diag_loss_write");

    file.append("repetition", repetition);
    file.append("s", s_ref_particle);
    file.append("s_n", sn_ref_particle);

    file.write_collective("coordinates_" + std::to_string(s_ref_particle), coords);
}

