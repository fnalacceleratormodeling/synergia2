#ifndef BUNCH_H_
#define BUNCH_H_

#include <mpi.h>
#include "utils/multi_array_typedefs.h"
#include "components/foundation/reference_particle.h"
#include "utils/commxx.h"

class Bunch;

class Fixed_t_z_converter
{
public:
    virtual void
    fixed_t_to_fixed_z(Bunch &bunch) = 0;
    virtual void
    fixed_z_to_fixed_t(Bunch &bunch) = 0;
};

class Fixed_t_z_zeroth : public Fixed_t_z_converter
{
public:
    void
    fixed_t_to_fixed_z(Bunch &bunch);
    void
    fixed_z_to_fixed_t(Bunch &bunch);
};

class Fixed_t_z_ballistic : public Fixed_t_z_converter
{
public:
    void
    fixed_t_to_fixed_z(Bunch &bunch);
    void
    fixed_z_to_fixed_t(Bunch &bunch);
};

class Bunch
{
public:
    enum State
    {
        fixed_z = 1, fixed_t = 2
    };
    static const int x = 0;
    static const int xp = 1;
    static const int y = 2;
    static const int yp = 3;
    static const int z = 4;
    static const int zp = 5;
    static const int t = 4;
    static const int tp = 5;
    static const int id = 6;
private:
    Reference_particle reference_particle;
    int particle_charge;
    MArray2d *local_particles;
    int local_num, total_num;
    double real_num;
    State state;
    bool particles_valid;
    Commxx comm;
    Fixed_t_z_converter *converter_ptr;
    Fixed_t_z_zeroth default_converter;
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
    void
    set_converter(Fixed_t_z_converter &converter);
    void
    convert_to_state(State state);
    Reference_particle &
    get_reference_particle();
    MArray2d_ref
    get_local_particles();
    int
    get_particle_charge();
    double
    get_mass();
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
