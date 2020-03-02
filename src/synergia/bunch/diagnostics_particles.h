#ifndef DIAGNOSTICS_PARTICLES_H_
#define DIAGNOSTICS_PARTICLES_H_

#include "synergia/bunch/diagnostics.h"

/// Diagnostics_particles dumps the state of particles in a bunch
class Diagnostics_particles : public Diagnostics
{

private:

    Bunch const* bunch_ptr;

    int num_part, offset;
    int num_spec_part, spec_offset;

public:

    Diagnostics_particles(
            int num_part = -1, int offset = 0,
            int num_spec_part = 0, int spec_offset = 0 )
        : Diagnostics("diagnostics_particles", false)
        , bunch_ptr(nullptr)
        , num_part(num_part) , offset(offset)
        , num_spec_part(num_spec_part), spec_offset(spec_offset)
    { }

private:

    void do_reduce(Commxx const& comm, int root) override { }
    void do_first_write(Hdf5_file& file) override { }

    void do_update(Bunch const& bunch) override 
    { bunch_ptr = &bunch; }

    void do_write(Hdf5_file& file) override
    { 
        assert(bunch_ptr != nullptr);
        bunch_ptr->write_file(file, num_part, offset, num_spec_part, spec_offset); 
        bunch_ptr = nullptr;
    }

    friend class cereal::access;

    template<class AR>
    void serialize(AR & ar)
    { 
        ar(cereal::base_class<Diagnostics>(this)); 
        ar(num_part);
        ar(offset);
        ar(num_spec_part);
        ar(spec_offset);
    }
};

CEREAL_REGISTER_TYPE(Diagnostics_particles);

#endif /* DIAGNOSTICS_PARTICLES_H_ */
