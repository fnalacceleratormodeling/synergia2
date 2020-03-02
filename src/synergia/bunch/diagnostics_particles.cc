#include "diagnostics_particles.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"

const char Diagnostics_particles::name[] = "diagnostics_particles";

Diagnostics_particles::Diagnostics_particles(
        std::string const& filename,
        int min_particle_id, 
        int max_particle_id, 
        std::string const& local_dir ) 
    : Diagnostics_particles::Diagnostics(Diagnostics_particles::name, filename, local_dir)
    , have_writers(false)
    , min_particle_id(min_particle_id)
    , max_particle_id(max_particle_id)
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
write_selected_particles(
        Hdf5_chunked_array2d_writer & writer,
        MArray2d_ref const & particles, 
        int local_num, 
        int min_particle_id,
        int max_particle_id )
{
    if (local_num == 0) 
    {
        return;
    }

    if ((min_particle_id == 0) && (max_particle_id == 0)) 
    {
        writer.write_chunk(
                particles[boost::indices[range(0, local_num)][range()]]);
    } 
    else 
    {
        for (int part = 0; part < local_num; ++part) 
        {
            int particle_id = int(particles[part][Bunch::id]);

            if ( (particle_id >= min_particle_id) && 
                 (particle_id <= max_particle_id)) 
            {
                writer.write_chunk(
                        particles[boost::indices[range(part, part + 1)][range()]]);
            }
        }
    }
}

void
Diagnostics_particles::receive_other_local_particles(
        std::vector<int> const & local_nums, 
        std::vector<int> const & local_nums_slots, 
        MArray2d_ref local_particles,
        int min_part_id,
        int max_part_id,
        Hdf5_file_sptr file_sptr)
{
    if (get_bunch().get_comm().has_this_rank())
    {
        int myrank = get_bunch().get_comm().get_rank();
        int size = get_bunch().get_comm().get_size();

        Hdf5_chunked_array2d_writer
                writer_particles(
                    file_sptr->get_h5file(),
                    "particles",
                    local_particles );

        for (int rank = 0; rank < size; ++rank) 
        {
            int local_num = local_nums[rank];
            int local_num_slots = local_nums_slots[rank];

            if (rank == myrank) 
            {
                write_selected_particles(
                        writer_particles, local_particles, 
                        local_num, min_part_id, max_part_id );
            } 
            else 
            {
                MPI_Status status;
                Raw_MArray2d received(
                        boost::extents[local_num_slots][7], 
                        boost::fortran_storage_order());

                int message_size = 7 * local_num_slots;
                MPI_Comm comm = get_bunch().get_comm().get();

                int error = MPI_Recv((void*) received.m.origin(), message_size,
                                     MPI_DOUBLE, rank, rank, comm, &status);

                if (error != MPI_SUCCESS) 
                {
                    throw std::runtime_error(
                                "Diagnostics_particles::receive_other_local_particles: MPI_Recv failed.");
                }

                write_selected_particles(
                        writer_particles, received.m, 
                        local_num, min_part_id, max_part_id );
            }
        }
    }
}

void
Diagnostics_particles::send_local_particles(
        int local_num_slots,
        MArray2d_ref local_particles,
        Diagnostics_write_helper & helper )
{
    if (get_bunch().get_comm().has_this_rank())
    {
#if 0
        int local_num_slots = get_bunch().get_local_num_slots();
        void * send_buffer = (void*) get_bunch().get_local_particles().origin();
#endif
        void * send_buffer = (void*)local_particles.origin();

        int message_size = 7 * local_num_slots;
        int receiver = helper.get_writer_rank();
        int rank = get_bunch().get_comm().get_rank();
        MPI_Comm comm = get_bunch().get_comm().get();

        int status = MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank, comm);

        if (status != MPI_SUCCESS) 
        {
            throw std::runtime_error(
                        "Diagnostics_particles::send_local_particles: MPI_Send failed.");
        }
    }
}

void
Diagnostics_particles::write()
{
    if (!get_bunch().get_comm().has_this_rank())
        return;

    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
    MPI_Comm comm = get_bunch().get_comm().get();
    int num_procs = get_bunch().get_comm().get_size();

    int status;

    {
        // local num and local num array size (slots)
        int local_num = get_bunch().get_local_num();
        int local_num_slots = get_bunch().get_local_num_slots();

        std::vector<int > local_nums(num_procs);
        std::vector<int > local_nums_slots(num_procs);

        void * local_nums_buf = (void *) &local_nums[0];
        void * local_nums_slots_buf = (void *) &local_nums_slots[0];

        int root = get_write_helper().get_writer_rank();

        status = MPI_Gather(
                (void*) &local_num, 1, MPI_INT, local_nums_buf, 
                1, MPI_INT, root, comm);

        if (status != MPI_SUCCESS) 
        {
            throw std::runtime_error(
                        "Diagnostics_particles::write: MPI_Gather local_nums failed.");
        }

        status = MPI_Gather(
                (void*) &local_num_slots, 1, MPI_INT, local_nums_slots_buf, 
                1, MPI_INT, root, comm);

        if (status != MPI_SUCCESS) 
        {
            throw std::runtime_error(
                        "Diagnostics_particles::write: MPI_Gather local_num_slots failed.");
        }

        // bunch particles
        if (get_write_helper().write_locally()) 
        {
            Hdf5_file_sptr file_sptr = get_write_helper().get_hdf5_file_sptr();

            receive_other_local_particles(
                    local_nums, 
                    local_nums_slots, 
                    get_bunch().get_local_particles(),
                    min_particle_id,
                    max_particle_id,
                    file_sptr );

            Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );

            int chg = get_bunch().get_reference_particle().get_charge();
            file_sptr->write(chg, "charge");

            double pmass = fourp.get_mass();
            file_sptr->write(pmass, "mass");

            double pz = fourp.get_momentum();
            file_sptr->write(pz, "pz");

            double tlen = get_bunch().get_reference_particle().get_s();
            file_sptr->write(tlen, "tlen");

            int rep = get_bunch().get_reference_particle().get_repetition();
            file_sptr->write(rep, "rep");

            double s_n = get_bunch().get_reference_particle().get_s_n();
            file_sptr->write(s_n, "s_n");
            file_sptr->write(0, "particles_storage_order");

            get_write_helper().finish_write();
        } 
        else 
        {
            send_local_particles(
                    local_num_slots, 
                    get_bunch().get_local_particles(), 
                    get_write_helper() );
        }
    }

    // any spectator particles?
    if (get_bunch().get_total_spectator_num() == 0)
        return;

    {
        // create/get write helper for the spectator particles
        Diagnostics_write_helper & helper = get_extra_write_helper("spectator");

        // local spectator num and local spectator array size (slots)
        int local_s_num = get_bunch().get_local_spectator_num();
        int local_s_num_slots = get_bunch().get_local_spectator_num_slots();

        std::vector<int> local_s_nums(num_procs);
        std::vector<int> local_s_nums_slots(num_procs);

        void * local_s_nums_buf = (void *) &local_s_nums[0];
        void * local_s_nums_slots_buf = (void *) &local_s_nums_slots[0];

        int s_root = helper.get_writer_rank();

        status = MPI_Gather(
                (void*) &local_s_num, 1, MPI_INT, local_s_nums_buf, 
                1, MPI_INT, s_root, comm);

        if (status != MPI_SUCCESS) 
        {
            throw std::runtime_error(
                        "Diagnostics_particles::write: MPI_Gather local_s_nums failed.");
        }

        status = MPI_Gather(
                (void*) &local_s_num_slots, 1, MPI_INT, local_s_nums_slots_buf, 
                1, MPI_INT, s_root, comm);

        if (status != MPI_SUCCESS) 
        {
            throw std::runtime_error(
                        "Diagnostics_particles::write: MPI_Gather local_s_num_slots failed.");
        }


        if (helper.write_locally()) 
        {
            Hdf5_file_sptr file_sptr = helper.get_hdf5_file_sptr();

            receive_other_local_particles(
                    local_s_nums, 
                    local_s_nums_slots, 
                    get_bunch().get_local_spectator_particles(),
                    0,  // always write out the entire block of spectator particles
                    0,
                    file_sptr );

            Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );

            int chg = get_bunch().get_reference_particle().get_charge();
            file_sptr->write(chg, "charge");

            double pmass = fourp.get_mass();
            file_sptr->write(pmass, "mass");

            double pz = fourp.get_momentum();
            file_sptr->write(pz, "pz");

            double tlen = get_bunch().get_reference_particle().get_s();
            file_sptr->write(tlen, "tlen");

            int rep = get_bunch().get_reference_particle().get_repetition();
            file_sptr->write(rep, "rep");

            double s_n = get_bunch().get_reference_particle().get_s_n();
            file_sptr->write(s_n, "s_n");
            file_sptr->write(0, "particles_storage_order");

            helper.finish_write();
        } 
        else 
        {
            send_local_particles(
                    local_s_num_slots,
                    get_bunch().get_local_spectator_particles(),
                    helper );
        }
    }
}

template<class Archive>
    void
    Diagnostics_particles::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
        ar & BOOST_SERIALIZATION_NVP(have_writers);
        ar & BOOST_SERIALIZATION_NVP(min_particle_id);
        ar & BOOST_SERIALIZATION_NVP(max_particle_id);
    }

template
void
Diagnostics_particles::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_particles::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_particles::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_particles::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_particles::~Diagnostics_particles()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_particles)
