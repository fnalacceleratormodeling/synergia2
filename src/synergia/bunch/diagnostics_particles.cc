#include "diagnostics_particles.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"

const char Diagnostics_particles::name[] = "diagnostics_particles";

Diagnostics_particles::Diagnostics_particles(std::string const& filename,
        int min_particle_id, int max_particle_id, std::string const& local_dir) :
        Diagnostics_particles::Diagnostics(Diagnostics_particles::name,
                filename, local_dir), have_writers(false), min_particle_id(
                min_particle_id), max_particle_id(max_particle_id)
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
    if (get_bunch().get_comm().has_this_rank()){
	int myrank = get_bunch().get_comm().get_rank();
	int size = get_bunch().get_comm().get_size();
	Hdf5_chunked_array2d_writer
		writer_particles(
			&(file_sptr->get_h5file()),
			"particles",
			get_bunch().get_local_particles());
	for (int rank = 0; rank < size; ++rank) {
	    int local_num = local_nums[rank];
	    if (rank == myrank) {
		write_selected_particles(writer_particles,
			get_bunch().get_local_particles(), local_num,
			min_particle_id, max_particle_id);
	    } else {
		MPI_Status status;
		Raw_MArray2d received(boost::extents[local_num][7]);
		int message_size = 7 * local_num;
		MPI_Comm comm = get_bunch().get_comm().get();
		int error = MPI_Recv((void*) received.m.origin(), message_size,
			MPI_DOUBLE, rank, rank, comm, &status);
		if (error != MPI_SUCCESS) {
		    throw std::runtime_error(
			    "Diagnostics_particles::receive_other_local_particles: MPI_Recv failed.");
		}
		write_selected_particles(writer_particles, received.m, local_num,
			min_particle_id, max_particle_id);
	    }
	}
    }
}

void
Diagnostics_particles::send_local_particles()
{
    if (get_bunch().get_comm().has_this_rank()){
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
}

void
Diagnostics_particles::write()
{
    if (get_bunch().get_comm().has_this_rank()){
      
	get_bunch().convert_to_state(get_bunch().fixed_z_lab);
	MPI_Comm comm = get_bunch().get_comm().get();
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
	    Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );
	    int chg = get_bunch().get_reference_particle().get_charge();
	    file_sptr->write(chg, "charge");
	    double pmass = fourp.get_mass();
	    file_sptr->write(pmass, "mass");
	    double pz = fourp.get_momentum();
	    file_sptr->write(pz, "pz");
	    double tlen =
		    get_bunch().get_reference_particle().get_s();
	    file_sptr->write(tlen, "tlen");
	    int rep = get_bunch().get_reference_particle().get_repetition();
	    file_sptr->write(rep, "rep");
	    double s_n = get_bunch().get_reference_particle().get_s_n();
	    file_sptr->write(s_n, "s_n");
	    get_write_helper().finish_write();
	} else {
	    send_local_particles();
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
