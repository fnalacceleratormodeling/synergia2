
#include "diagnostics_particles.h"
#include "synergia/bunch/bunch.h"

Diagnostics_particles::Diagnostics_particles(
        int num_part, 
        int offset,
        int num_spec_part, 
        int spec_offset)
    : Diagnostics("diagnostics_particles", false)
    , bunch_ptr(nullptr)
    , num_part(num_part) , offset(offset)
    , num_spec_part(num_spec_part), spec_offset(spec_offset)
{ }


void Diagnostics_particles::do_write(Hdf5_file& file)
{ 
    assert(bunch_ptr != nullptr);
    bunch_ptr->write_file(file, num_part, offset, num_spec_part, spec_offset); 
    bunch_ptr = nullptr;
}

