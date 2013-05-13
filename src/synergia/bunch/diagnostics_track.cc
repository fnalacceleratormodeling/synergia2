#include "diagnostics_track.h"


const char Diagnostics_track::name[] = "diagnostics_track";

Diagnostics_track::Diagnostics_track(std::string const& filename,
        int particle_id, std::string const& local_dir) :
        Diagnostics(Diagnostics_track::name, filename, local_dir), have_writers(
                false), found(false), first_search(true), last_index(-1), particle_id(
                particle_id), coords(boost::extents[6])
{
}

Diagnostics_track::Diagnostics_track() : have_writers(false)
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
            get_bunch().get_comm_sptr(), get_local_dir(),
            get_bunch().get_comm().get_rank());
}

void
Diagnostics_track::update()
{
    if (get_bunch().get_comm().has_this_rank()){
	get_bunch().convert_to_state(get_bunch().fixed_z_lab);
	repetition = get_bunch().get_reference_particle().get_repetition();
	trajectory_length
		= get_bunch().get_reference_particle().get_trajectory_length();
	if (found || first_search) {
	    int index = 0;
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
		s_n = get_bunch().get_reference_particle().get_s_n();
		repetition = get_bunch().get_reference_particle().get_repetition();
		trajectory_length
			= get_bunch().get_reference_particle().get_trajectory_length();
	    }
	    first_search = false;
	}
    }
}

void
Diagnostics_track::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );
        int chg = get_bunch().get_reference_particle().get_charge();
        file_sptr->write(chg, "charge");
        double pmass = fourp.get_mass();
        file_sptr->write(pmass, "mass");
        double pz = fourp.get_momentum();
        file_sptr->write(pz, "pz");

        writer_coords = new Hdf5_serial_writer<MArray1d_ref > (file_sptr,
                "coords");
        writer_s_n = new Hdf5_serial_writer<double > (file_sptr, "s_n");
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
    if (get_bunch().get_comm().has_this_rank()){
	get_bunch().convert_to_state(get_bunch().fixed_z_lab);
	if (found) {
	    init_writers(get_write_helper().get_hdf5_file_sptr());
	    writer_coords->append(coords);
	    writer_s_n->append(s_n);
	    writer_repetition->append(repetition);
	    writer_trajectory_length->append(trajectory_length);
	    get_write_helper().finish_write();
	}
    }
}

template<class Archive>
    void
    Diagnostics_track::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                & BOOST_SERIALIZATION_NVP(have_writers)
                & BOOST_SERIALIZATION_NVP(found)
                & BOOST_SERIALIZATION_NVP(first_search)
                & BOOST_SERIALIZATION_NVP(last_index)
                & BOOST_SERIALIZATION_NVP(particle_id)
                & BOOST_SERIALIZATION_NVP(s_n)
                & BOOST_SERIALIZATION_NVP(writer_s_n)
                & BOOST_SERIALIZATION_NVP(repetition)
                & BOOST_SERIALIZATION_NVP(writer_repetition)
                & BOOST_SERIALIZATION_NVP(trajectory_length)
                & BOOST_SERIALIZATION_NVP(writer_trajectory_length)
                & BOOST_SERIALIZATION_NVP(coords)
                & BOOST_SERIALIZATION_NVP(writer_coords);
    }

template
void
Diagnostics_track::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_track::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_track::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_track::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_track::~Diagnostics_track()
{
    if (have_writers) {
        delete writer_trajectory_length;
        delete writer_repetition;
        delete writer_s_n;
        delete writer_coords;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_track)

