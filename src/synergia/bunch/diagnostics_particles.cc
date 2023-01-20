
#include "diagnostics_particles.h"
#include "synergia/bunch/bunch.h"
#include "synergia_config.h"

Diagnostics_particles::Diagnostics_particles(std::string const& filename,
                                             int num_part,
                                             int offset,
                                             int num_spec_part,
                                             int spec_offset)
#ifdef SYNERGIA_HAVE_OPENPMD
    // use a single file and iterations with openPMD
    : Diagnostics("diagnostics_particles", filename, true)
#else
    // use a new file for each write with old HDF5 backend
    : Diagnostics("diagnostics_particles", filename, false)
#endif
    , bunch_ref(std::nullopt)
    , num_part(num_part)
    , offset(offset)
    , num_spec_part(num_spec_part)
    , spec_offset(spec_offset)
{}

void
Diagnostics_particles::do_first_write(io_device& file)
{

#ifdef SYNERGIA_HAVE_OPENPMD
    assert(bunch_ref.has_value());
    auto const& ref_part = bunch_ref.value().get().get_reference_particle();
    file.setAttribute("charge", ref_part.get_charge());
    file.setAttribute("mass", ref_part.get_mass());
    file.setAttribute("s", ref_part.get_s());
    file.setAttribute("s_n", ref_part.get_s_n());
    file.setAttribute("repetition", ref_part.get_repetition());
    file.setAttribute("pz", ref_part.get_momentum());
#else
    // nothing to do if using the old HDF5 backend!

#endif
    return;
}

void
Diagnostics_particles::do_write(io_device& file, const size_t iteration)
{

#ifdef SYNERGIA_HAVE_OPENPMD
    auto i = file.iterations[iteration];

#else
    assert(bunch_ref.has_value());
    auto const& ref_part = bunch_ref.value().get().get_reference_particle();

    file.write("charge", ref_part.get_charge());
    file.write("mass", ref_part.get_four_momentum().get_mass());

    file.write("s", ref_part.get_s());
    file.write("s_n", ref_part.get_s_n());
    file.write("repetition", ref_part.get_repetition());
    file.write("pz", ref_part.get_momentum());

    bunch_ref.value().get().write_file(
        file, num_part, offset, num_spec_part, spec_offset);

    // reset bunch_ref
    bunch_ref = std::nullopt;
#endif
}
