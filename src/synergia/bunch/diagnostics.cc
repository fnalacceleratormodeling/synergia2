#include "diagnostics.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "synergia/utils/eigen2/Eigen/Core"
#include "synergia/utils/eigen2/Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

Diagnostics::Diagnostics(std::string const& name, std::string const& filename) :
    name(name), filename(filename), have_bunch_(false), write_helper_ptr(0),
            have_write_helper_(0)
{
}

std::string const&
Diagnostics::get_filename() const
{
    return filename;
}

void
Diagnostics::set_bunch_sptr(Bunch_sptr bunch_sptr)
{
    this->bunch_sptr = bunch_sptr;
    have_bunch_ = true;
}

bool
Diagnostics::have_bunch() const
{
    return have_bunch_;
}

Bunch &
Diagnostics::get_bunch()
{
    if (!have_bunch_) {
        throw std::runtime_error(name + ": bunch not set");
    }
    return *bunch_sptr;
}

void
Diagnostics::delete_write_helper_ptr()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
        have_write_helper_ = false;
    }
}

Diagnostics_write_helper *
Diagnostics::new_write_helper_ptr()
{
    delete_write_helper_ptr();
    return new Diagnostics_write_helper(get_filename(),
            is_serial(), get_bunch().get_comm());
}

bool
Diagnostics::have_write_helper() const
{
    return have_write_helper_;
}

Diagnostics_write_helper &
Diagnostics::get_write_helper()
{
    if (!have_write_helper_) {
        write_helper_ptr = new_write_helper_ptr();
        have_write_helper_ = true;
    }
    return *write_helper_ptr;
}

Diagnostics::Diagnostics()
{
}

Diagnostics::~Diagnostics()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
    }
}

BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics)



const char Diagnostics_particles::name[] = "diagnostics_particles";

Diagnostics_particles::Diagnostics_particles(
        std::string const& filename, int min_particle_id, int max_particle_id,
        int write_skip) :
    Diagnostics_particles::Diagnostics(Diagnostics_particles::name, filename),
            min_particle_id(min_particle_id), max_particle_id(max_particle_id),
            have_writers(false)
{
}

Diagnostics_particles::Diagnostics_particles()
{
}

bool
Diagnostics_particles::is_serial() const
{
    return false;
}

void
Diagnostics_particles::update()
{
}

// write_selected_particles is a local function
void
write_selected_particles(Hdf5_chunked_array2d_writer & writer,
        MArray2d_ref const & particles, int local_num, int min_particle_id,
        int max_particle_id)
{

    if ((min_particle_id == 0) && (max_particle_id == 0)) {
        writer.write_chunk(
                particles[boost::indices[range(0, local_num)][range()]]);
    } else {
        for (int part = 0; part < local_num; ++part) {
            int particle_id = int(particles[part][Bunch::id]);
            if ((particle_id >= min_particle_id) && (particle_id
                    <= max_particle_id)) {
                writer.write_chunk(
                        particles[boost::indices[range(part, part + 1)][range()]]);
            }
        }
    }
}

void
Diagnostics_particles::receive_other_local_particles(
        std::vector<int > const& local_nums, Hdf5_file_sptr file_sptr)
{
    int myrank = get_bunch().get_comm().get_rank();
    int size = get_bunch().get_comm().get_size();
    Hdf5_chunked_array2d_writer
            writer_particles(
                    &(file_sptr->get_h5file()),
                    "particles",
                    get_bunch().get_local_particles()[boost::indices[range(0, 1)][range()]]);
    for (int rank = 0; rank < size; ++rank) {
        int local_num = local_nums[rank];
        if (rank == myrank) {
            write_selected_particles(writer_particles,
                    get_bunch().get_local_particles(), local_num,
                    min_particle_id, max_particle_id);
        } else {
            MPI_Status status;
            MArray2d received(boost::extents[local_num][7]);
            int message_size = 7 * local_num;
            MPI_Comm comm = get_bunch().get_comm().get();
            MPI_Recv((void*) received.origin(), message_size, MPI_DOUBLE, rank,
                    rank, comm, &status);
            if (status.MPI_ERROR != MPI_SUCCESS) {
                throw std::runtime_error(
                        "Diagnostics_particles::receive_other_local_particles: MPI_Recv failed.");
            }
            write_selected_particles(writer_particles, received, local_num,
                    min_particle_id, max_particle_id);
        }
    }
}

void
Diagnostics_particles::send_local_particles()
{
    int local_num = get_bunch().get_local_num();
    void * send_buffer =
            (void*) get_bunch().get_local_particles()[boost::indices[range(0,
                    local_num)][range()]].origin();
    int status;
    int message_size = 7 * local_num;
    int receiver = get_write_helper().get_writer_rank();
    int rank = get_bunch().get_comm().get_rank();
    MPI_Comm comm = get_bunch().get_comm().get();
    status = MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
            comm);
    if (status != MPI_SUCCESS) {
        throw std::runtime_error(
                "Diagnostics_particles::send_local_particles: MPI_Send failed.");
    }
}

void
Diagnostics_particles::write()
{
    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
    int writer_rank = get_write_helper().get_writer_rank();
    MPI_Comm comm = get_bunch().get_comm().get();
    int icount;
    icount = get_write_helper().get_count();
    MPI_Bcast((void *) &icount, 1, MPI_INT, writer_rank, comm);

    //    std::cout<<" icount ="<<icount<<"  count="<< get_write_helper().get_count() <<" on rank ="<< get_bunch().get_comm().get_rank()<<std::endl;

    if (icount % get_write_helper().get_iwrite_skip() != 0) {
        if (get_write_helper().write_locally()) get_write_helper().increment_count();
        return;
    }

    int local_num = get_bunch().get_local_num();
    int num_procs = get_bunch().get_comm().get_size();
    std::vector<int > local_nums(num_procs);
    void * local_nums_buf = (void *) &local_nums[0];
    int root = get_write_helper().get_writer_rank();
    int status;
    status = MPI_Gather((void*) &local_num, 1, MPI_INT, local_nums_buf, 1,
            MPI_INT, root, comm);
    if (status != MPI_SUCCESS) {
        throw std::runtime_error(
                "Diagnostics_particles::write: MPI_Gather failed.");
    }
    if (get_write_helper().write_locally()) {
        Hdf5_file_sptr file_sptr = get_write_helper().get_hdf5_file_sptr();
        receive_other_local_particles(local_nums, file_sptr);
        double pz = get_bunch().get_reference_particle().get_momentum();
        file_sptr->write(pz, "pz");
        double tlen =
                get_bunch().get_reference_particle().get_trajectory_length();
        file_sptr->write(tlen, "tlen");
        int rep = get_bunch().get_reference_particle().get_repetition();
        file_sptr->write(rep, "rep");
        double s = get_bunch().get_reference_particle().get_s();
        file_sptr->write(s, "s");
        get_write_helper().finish_write();
    } else {
        send_local_particles();
    }
}

Diagnostics_particles::~Diagnostics_particles()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_particles)

const char Diagnostics_track::name[] = "diagnostics_track";

Diagnostics_track::Diagnostics_track(
        std::string const& filename, int particle_id) :
    Diagnostics(Diagnostics_track::name, filename), have_writers(false),
            coords(boost::extents[6]), found(false), first_search(true),
            particle_id(particle_id), last_index(-1)
{
}

Diagnostics_track::Diagnostics_track()
{
}

bool
Diagnostics_track::is_serial() const
{
    return true;
}

Diagnostics_write_helper *
Diagnostics_track::new_write_helper_ptr()
{
    delete_write_helper_ptr();
    return new Diagnostics_write_helper(get_filename(), true,
            get_bunch().get_comm(), get_bunch().get_comm().get_rank());
}

void
Diagnostics_track::update()
{
    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
    repetition = get_bunch().get_reference_particle().get_repetition();
    trajectory_length
            = get_bunch().get_reference_particle().get_trajectory_length();
    if (found || first_search) {
        int index;
        found = false;
        if ((last_index > -1) && (last_index < get_bunch().get_local_num())) {
            if (particle_id
                    == static_cast<int > (get_bunch().get_local_particles()[Bunch::id][last_index])) {
                index = last_index;
                found = true;
            }
        }
        if (!found) {
            index = 0;
            while ((index < get_bunch().get_local_num())
                    && (particle_id
                            != static_cast<int > (get_bunch().get_local_particles()[index][Bunch::id]))) {
                index += 1;
            }
            if (index < get_bunch().get_local_num()) {
                found = true;
            } else {
                found = false;
            }
        }
        if (found) {
            if (first_search) {
                get_write_helper();
            }
            coords[0] = get_bunch().get_local_particles()[index][0];
            coords[1] = get_bunch().get_local_particles()[index][1];
            coords[2] = get_bunch().get_local_particles()[index][2];
            coords[3] = get_bunch().get_local_particles()[index][3];
            coords[4] = get_bunch().get_local_particles()[index][4];
            coords[5] = get_bunch().get_local_particles()[index][5];
            s = get_bunch().get_reference_particle().get_s();
            repetition = get_bunch().get_reference_particle().get_repetition();
            trajectory_length
                    = get_bunch().get_reference_particle().get_trajectory_length();
        }
        first_search = false;
    }
}

void
Diagnostics_track::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_coords = new Hdf5_serial_writer<MArray1d_ref > (file_sptr,
                "coords");
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");
        writer_repetition = new Hdf5_serial_writer<int > (file_sptr,
                "repetition");
        writer_trajectory_length = new Hdf5_serial_writer<double > (file_sptr,
                "trajectory_length");
        have_writers = true;
    }
}

void
Diagnostics_track::write()
{
    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
    if (found) {
        init_writers(get_write_helper().get_hdf5_file_sptr());
        writer_coords->append(coords);
        writer_s->append(s);
        writer_repetition->append(repetition);
        writer_trajectory_length->append(trajectory_length);
        get_write_helper().finish_write();
    }
}

Diagnostics_track::~Diagnostics_track()
{
    if (have_writers) {
        delete writer_trajectory_length;
        delete writer_repetition;
        delete writer_s;
        delete writer_coords;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_track)

const char Diagnostics_reference_particle::name[] = "diagnostics_reference_particle";

Diagnostics_reference_particle::Diagnostics_reference_particle(
        std::string const& filename) :
            Diagnostics_reference_particle::Diagnostics(
                    Diagnostics_reference_particle::name, filename),
            have_writers(false), writer_beta(0), writer_gamma(0),
            writer_state(0), writer_s(0)
{
}

Diagnostics_reference_particle::Diagnostics_reference_particle()
{
}

bool
Diagnostics_reference_particle::is_serial() const
{
    return true;
}

void
Diagnostics_reference_particle::update()
{
}

void
Diagnostics_reference_particle::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_beta = new Hdf5_serial_writer<double > (file_sptr, "beta");
        writer_gamma = new Hdf5_serial_writer<double > (file_sptr, "gamma");
        writer_state = new Hdf5_serial_writer<MArray1d_ref > (file_sptr,
                "state");
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");
        have_writers = true;
    }
}

void
Diagnostics_reference_particle::write()
{
    if (get_write_helper().write_locally()) {
        init_writers(get_write_helper().get_hdf5_file_sptr());
        double beta = get_bunch().get_reference_particle().get_beta();
        writer_beta->append(beta);
        double gamma = get_bunch().get_reference_particle().get_gamma();
        writer_gamma->append(gamma);
        MArray1d state(get_bunch().get_reference_particle().get_state());
        writer_state->append(state);
        double s = get_bunch().get_reference_particle().get_s();
        writer_s->append(s);
        get_write_helper().finish_write();
    }
}

Diagnostics_reference_particle::~Diagnostics_reference_particle()
{
    if (have_writers) {
        delete writer_beta;
        delete writer_gamma;
        delete writer_state;
        delete writer_s;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_reference_particle)

