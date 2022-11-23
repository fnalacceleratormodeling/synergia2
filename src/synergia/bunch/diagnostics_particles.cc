
#include "diagnostics_particles.h"
#include "synergia/bunch/bunch.h"

Diagnostics_particles::Diagnostics_particles(std::string const& filename,
                                             int num_part,
                                             int offset,
                                             int num_spec_part,
                                             int spec_offset)
    : Diagnostics("diagnostics_particles", filename, false)
    , bunch_ptr(nullptr)
    , num_part(num_part)
    , offset(offset)
    , num_spec_part(num_spec_part)
    , spec_offset(spec_offset)
{}

void
Diagnostics_particles::do_write(Hdf5_file& file)
{
    assert(bunch_ptr != nullptr);
    auto const& ref = bunch_ptr->get_reference_particle();

    file.write("charge", ref.get_charge());
    file.write("mass", ref.get_four_momentum().get_mass());

    file.write("s", ref.get_s());
    file.write("s_n", ref.get_s_n());
    file.write("repetition", ref.get_repetition());
    file.write("pz", ref.get_momentum());

    bunch_ptr->write_file(file, num_part, offset, num_spec_part, spec_offset);
    bunch_ptr = nullptr;
}
