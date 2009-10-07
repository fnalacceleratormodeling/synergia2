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
    local_particles = new MArray2d(boost::extents[7][local_num]);
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
}
//int old_total_num = mbs.total_num;
// rank = 0;
// MPI_Comm_size(MPI_COMM_WORLD, &size);
// if (size > 1) {
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Allreduce(&mbs.local_num, &mbs.total_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
// } else {
//     mbs.total_num = mbs.local_num;
// }

void
Bunch::update_total_num()
{
    int old_total_num = total_num;
    MPI_Allreduce(&local_num, &total_num,1,MPI_INT,MPI_SUM,comm.get());
    real_num = (total_num*real_num)/old_total_num;
}
//	void set_total_num(int total_num, bool update_current=true);
//	void set_total_current(double total_current);
//
//	Array_2d<double>& get_local_particles();
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

//	double get_total_current();
//	Bunch_state get_state();
//
//	int convert_to_state(Bunch_state state);
//	void update_total_num();
//
//	void inject(Bunch bunch);

Bunch::~Bunch()
{
    delete local_particles;
}
