#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

#include "utils/multi_array_typedefs.h"
#include "components/foundation/four_momentum.h"

class Reference_particle
{
private:
    Four_momentum four_momentum;
    MArray1d state;
public:
    Reference_particle(Four_momentum const& four_momentum);
    Reference_particle(Four_momentum const& four_momentum,
            Const_MArray1d_ref state);

    void
    set_four_momentum(Four_momentum const& four_momentum);
    void
    set_state(Const_MArray1d_ref state);

    Four_momentum &
    get_four_momentum();
    Const_MArray1d_ref
    get_state();
};

#endif /* REFERENCE_PARTICLE_H_ */
