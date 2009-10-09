#ifndef BUNCH_H_
#define BUNCH_H_

#include <mpi.h>
#include "utils/multi_array_typedefs.h"
#include "components/foundation/reference_particle.h"
#include "utils/commxx.h"

class Bunch
{
public:
    enum State
    {
        fixed_z = 1, fixed_t = 2
    };

private:
    Reference_particle reference_particle;
    int particle_charge;
    MArray2d *local_particles;
    int local_num, total_num;
    double real_num;
    State state;
    bool particles_valid;
    Commxx comm;
public:
    Bunch(Reference_particle const& reference_particle, int particle_charge,
            int total_num, double real_num, Commxx const& comm);
    void
    set_particle_charge(int particle_charge);
    void
    set_real_num(double real_num);
    void
    set_local_num(int local_num);
    void
    update_total_num();
//    void
//    set_state(Bunch_state state);
    MArray2d_ref
    get_local_particles();
    int
    get_particle_charge();
    double
    get_real_num();
    int
    get_local_num();
    int
    get_total_num();
    State
    get_state();
    //	void inject(Bunch bunch);
    ~Bunch();
};

#endif /* BUNCH_H_ */
