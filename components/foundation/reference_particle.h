#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

#include "utils/multi_array_typedefs.h"

class Reference_particle
{
private:
    double total_energy;
    MArray1d state;
public:
    Reference_particle(double total_energy);
    Reference_particle(double total_energy, Const_MArray1d_ref state);

    void
    set_total_energy(double total_energy);
    void
    set_state(Const_MArray1d_ref state);

    double
    get_total_energy();
    MArray1d_ref
    get_state();
};

#endif /* REFERENCE_PARTICLE_H_ */
