#include "synergia/bunch/diagnostics_loss.h"
#include "synergia/bunch/bunch.h"
#include "synergia/utils/simple_timer.h"

void
Diagnostics_loss::do_update(Bunch const& bunch)
{
    scoped_simple_timer timer("diag_loss_update");

    auto ref = bunch.get_reference_particle();
    repetition = ref.get_repetition();
    s_ref_particle = ref.get_s();
    sn_ref_particle = ref.get_s_n();
    bucket_index = bunch.get_bucket_index();

    coords = bunch.get_particles_last_discarded();
}

void
Diagnostics_loss::do_first_write(io_device& file)
{
#ifdef SYNERGIA_HAVE_OPENPMD
    file.setAttribute("bucket_index", bucket_index);
#else
    file.append("bucket_index", bucket_index);
#endif
    return;
}

void
Diagnostics_loss::do_write(io_device& file, size_t iteration)
{
    scoped_simple_timer timer("diag_loss_write");
#ifdef SYNERGIA_HAVE_OPENPMD
    file.setAttribute("repetition", repetition);
    file.setAttribute("s", s_ref_particle);
    file.setAttribute("s_n", sn_ref_particle);
    // TODO: Write coords!

#else
    file.append("repetition", repetition);
    file.append("s", s_ref_particle);
    file.append("s_n", sn_ref_particle);

    file.write_collective("coordinates_" + std::to_string(s_ref_particle),
                          coords);
#endif

    return;
}
