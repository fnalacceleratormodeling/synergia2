#include "bunch.h"
#include "utils/parallel_utils.h"
#include <stdexcept>
//        double gamma = -ref_particle(5);
//        for (int i = 0; i < local_num; ++i) {
//            local_particles(0, i) /= units(0);
//            local_particles(2, i) /= units(2);
//            double xp = local_particles(1, i);
//            double yp = local_particles(3, i);
//            double rcp_gammai = 1.0 / (gamma - local_particles(5, i));
//            double betai = sqrt(1.0 - rcp_gammai * rcp_gammai * (1 + xp * xp + yp * yp));
//            local_particles(4, i) *= -gamma * betai / units(0); // units(0)

void
Fixed_t_z_zeroth::fixed_t_to_fixed_z(Bunch &bunch)
{
    std::cout << "zeroth fixed_t_to_fixed_z\n";
}

void
Fixed_t_z_zeroth::fixed_z_to_fixed_t(Bunch &bunch)
{
    //double gamma_ref = bunch.get_reference_particle().get_gamma();
}

void
Fixed_t_z_ballistic::fixed_t_to_fixed_z(Bunch &bunch)
{
    std::cout << "ballistic fixed_t_to_fixed_z\n";
}

void
Fixed_t_z_ballistic::fixed_z_to_fixed_t(Bunch &bunch)
{
    std::cout << "ballistic fixed_z_to_fixed_t\n";
}

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
    reference_particle(reference_particle), comm(comm), default_converter()
{
    this->particle_charge = particle_charge;
    this->total_num = total_num;
    this->real_num = real_num;
    local_num = decompose_1d_local(comm, total_num);
    local_particles = new MArray2d(boost::extents[local_num][7]);
    state = fixed_z;
    particles_valid = false;
    converter_ptr = &default_converter;
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

void
Bunch::set_converter(Fixed_t_z_converter &converter)
{
    this->converter_ptr = &converter;
}

void
Bunch::convert_to_state(State state)
{
    if (this->state != state) {
        if (state == fixed_t) {
            converter_ptr->fixed_z_to_fixed_t(*this);
            this->state = fixed_t;
        } else if (state == fixed_z) {
            converter_ptr->fixed_t_to_fixed_z(*this);
            this->state = fixed_z;
        } else {
            throw std::runtime_error("Unknown state in Bunch::convert_to_state");
        }
    }
}

Reference_particle &
Bunch::get_reference_particle()
{
    return reference_particle;
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
