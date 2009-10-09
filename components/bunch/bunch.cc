#include "bunch.h"
#include "utils/parallel_utils.h"
//
//    Reference_particle reference_particle;
//    int particle_charge;
//    MArray2d local_particles;
//    int local_num, total_num;
//    double real_num;
//    Bunch_state state;
//    bool particles_valid;
//    bool total_num_valid;
//    MPI_comm comm;


Bunch::Bunch(Reference_particle const& reference_particle, int particle_charge,
        int total_num, double real_num, Commxx const& comm) :
    reference_particle(reference_particle), comm(comm)
{
    this->particle_charge = particle_charge;
    this->total_num = total_num;
    this->real_num = real_num;
    local_num = decompose_1d_local(comm, total_num);
    local_particles = new MArray2d(boost::extents[local_num][7]);
    state = fixed_z;
    particles_valid = false;
}

void
Bunch::set_particle_charge(int particle_charge)
{
    this->particle_charge = particle_charge;
}

void
Bunch::set_real_num(double real_num)
{
    this->real_num = real_num;
}

void
Bunch::set_local_num(int local_num)
{
    this->local_num = local_num;
    if (local_particles->shape()[0] < local_num) {
        local_particles->resize(boost::extents[local_num][7]);
    }
}

void
Bunch::update_total_num()
{
    int old_total_num = total_num;
    MPI_Allreduce(&local_num, &total_num, 1, MPI_INT, MPI_SUM, comm.get());
    real_num = (total_num * real_num) / old_total_num;
}

MArray2d_ref
Bunch::get_local_particles()
{
    return *local_particles;
}

int
Bunch::get_particle_charge()
{
    return particle_charge;
}

double
Bunch::get_real_num()
{
    return real_num;
}

int
Bunch::get_local_num()
{
    return local_num;
}

int
Bunch::get_total_num()
{
    return total_num;
}

Bunch::State
Bunch::get_state()
{
    return state;
}

Bunch::~Bunch()
{
    delete local_particles;
}
