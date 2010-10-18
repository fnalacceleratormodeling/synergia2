#include "bunch.h"
#include "synergia/utils/parallel_utils.h"
#include <stdexcept>
#include <cmath>

const int Bunch::x;
const int Bunch::xp;
const int Bunch::y;
const int Bunch::yp;
const int Bunch::z;
const int Bunch::zp;
const int Bunch::t;
const int Bunch::tp;
const int Bunch::id;

class Particle_id_offset
{
private:
    int offset;
public:
    Particle_id_offset() :
        offset(0)
    {
    }

    int
    get(int request_num)
    {
        MPI_Bcast((void *) &offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int old_offset = offset;
        int total_num;
        MPI_Reduce((void*) &request_num, (void*) &total_num, 1, MPI_INT,
                MPI_SUM, 0, MPI_COMM_WORLD);
        offset += total_num;
        return old_offset;
    }

};

static Particle_id_offset particle_id_offset;

void
Bunch::assign_ids(int local_offset)
{
    int global_offset, request_num;
    if (comm.get_rank() == 0) {
        request_num = total_num;
    } else {
        request_num = 0;
    }
    global_offset = particle_id_offset.get(request_num);
    for (int i = 0; i < local_num; ++i) {
        (*local_particles)[i][id] = i + local_offset + global_offset;
    }
}

void
Bunch::construct(int particle_charge, int total_num, double real_num)
{
    this->particle_charge = particle_charge;
    this->total_num = total_num;
    this->real_num = real_num;
    state = fixed_z;
    converter_ptr = &default_converter;
    local_num = decompose_1d_local(comm, total_num);
    std::vector<int > offsets(comm.get_size()), counts(comm.get_size());
    decompose_1d(comm, total_num, offsets, counts);
    local_num = counts[comm.get_rank()];
    local_particles = new MArray2d(boost::extents[local_num][7]);
    assign_ids(offsets[comm.get_rank()]);
}

Bunch::Bunch(Reference_particle const& reference_particle, int total_num,
        double real_num, Commxx const& comm) :
    reference_particle(reference_particle), comm(comm), default_converter()
{
    construct(reference_particle.get_charge(), total_num, real_num);
}

Bunch::Bunch(Reference_particle const& reference_particle, int total_num,
        double real_num, Commxx const& comm, int particle_charge) :
    reference_particle(reference_particle), comm(comm), default_converter()
{
    construct(particle_charge, total_num, real_num);
}

Bunch::Bunch(Bunch const& bunch) :
    reference_particle(bunch.reference_particle), comm(bunch.comm),
            default_converter()
{
    particle_charge = bunch.particle_charge;
    total_num = bunch.total_num;
    real_num = bunch.real_num;
    local_num = bunch.local_num;
    local_particles = new MArray2d(*(bunch.local_particles));
    state = bunch.state;
    if (bunch.converter_ptr == &(bunch.default_converter)) {
        converter_ptr = &default_converter;
    } else {
        converter_ptr = bunch.converter_ptr;
    }
}

Bunch &
Bunch::operator=(Bunch const& bunch)
{
    if (this != &bunch) {
        reference_particle = bunch.reference_particle;
        comm = bunch.comm;
        particle_charge = bunch.particle_charge;
        total_num = bunch.total_num;
        real_num = bunch.real_num;
        local_num = bunch.local_num;
        local_particles = new MArray2d(*(bunch.local_particles));
        state = bunch.state;
        if (bunch.converter_ptr == &(bunch.default_converter)) {
            converter_ptr = &default_converter;
        } else {
            converter_ptr = bunch.converter_ptr;
        }
    }
    return *this;
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
    if (local_num > this->local_num) {
        throw std::runtime_error(
                "set_local_num can only decrease the number of particles");
    }
    this->local_num = local_num;
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

Reference_particle const&
Bunch::get_reference_particle() const
{
    return reference_particle;
}

MArray2d_ref
Bunch::get_local_particles()
{
    return *local_particles;
}

Const_MArray2d_ref
Bunch::get_local_particles() const
{
    return *local_particles;
}

int
Bunch::get_particle_charge() const
{
    return particle_charge;
}

double
Bunch::get_mass() const
{
    return reference_particle.get_four_momentum().get_mass();
}

double
Bunch::get_real_num() const
{
    return real_num;
}

int
Bunch::get_local_num() const
{
    return local_num;
}

int
Bunch::get_total_num() const
{
    return total_num;
}

Bunch::State
Bunch::get_state() const
{
    return state;
}

Commxx const&
Bunch::get_comm() const
{
    return comm;
}

void
Bunch::inject(Bunch const& bunch)
{
    const double weight_tolerance = 1.0e-10;
    if (std::abs(real_num / total_num - bunch.get_real_num()
            / bunch.get_total_num()) > weight_tolerance) {
        throw std::runtime_error(
                "Bunch.inject: macroparticle weight of injected bunch does not match.");
    }
    int old_local_num = local_num;
    local_num += bunch.get_local_num();
    Const_MArray2d_ref injected_particles(bunch.get_local_particles());
    MArray1d ref_state_diff(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        ref_state_diff[i] = bunch.get_reference_particle().get_state()[i]
                - reference_particle.get_state()[i];
    }
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            (*local_particles)[old_local_num + part][i]
                    = injected_particles[part][i] + ref_state_diff[i];
        }
        (*local_particles)[old_local_num + part][Bunch::id]
                = injected_particles[part][Bunch::id];
    }
    update_total_num();
}

Bunch::~Bunch()
{
    delete local_particles;
}
