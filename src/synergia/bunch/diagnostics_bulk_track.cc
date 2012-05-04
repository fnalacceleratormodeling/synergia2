#include <string>
#include <iostream>
#include "diagnostics_bulk_track.h"

const char Diagnostics_bulk_track::name[] = "diagnostics_bulk_track";

Diagnostics_bulk_track::Diagnostics_bulk_track(
    std::string const& filename, int num_tracks, bool iocc_verbose) :
  Diagnostics(Diagnostics_bulk_track::name, filename), num_tracks(num_tracks),
  iocc_verbose(iocc_verbose),
  have_writers(false), first_search(true), diag_track_status(),
  track_coords(boost::extents[num_tracks][7])
  
{
}

Diagnostics_bulk_track::Diagnostics_bulk_track()
{
}

bool
Diagnostics_bulk_track::is_serial() const
{
    return true;
}

Diagnostics_write_helper *
Diagnostics_bulk_track::new_write_helper_ptr()
{
    delete_write_helper_ptr();
    std::stringstream sstream;
    std::string filename(get_filename());
    int idx = filename.rfind('.');
    std::string filename_base, filename_suffix;
    if (idx == std::string::npos) {
      filename_base = filename;
      filename_suffix = "";
    } else {
      filename_base = filename.substr(0,idx);
      filename_suffix = filename.substr(idx);
    }

    sstream << filename_base;
    sstream << "_";
    sstream << std::setw(4);

    sstream << std::setfill('0');
    sstream << get_bunch().get_comm().get_rank();
    sstream << filename_suffix;
    return new Diagnostics_write_helper(sstream.str(), true,
            get_bunch().get_comm(), get_bunch().get_comm().get_rank());
}

void
Diagnostics_bulk_track::update()
{
  if (diag_track_status.empty()) {
    if (num_tracks > get_bunch().get_local_num()) {
      num_tracks = get_bunch().get_local_num();
    }
    for (int idxtrk=0; idxtrk<num_tracks; ++idxtrk) {
      Track_status dts;
      dts.found = true;
      dts.last_index = idxtrk;
      dts.particle_id = static_cast<int > (get_bunch().get_local_particles()[idxtrk][Bunch::id]);
      diag_track_status.push_back(dts);
    }
  }
  get_bunch().convert_to_state(get_bunch().fixed_z_lab);
  repetition = get_bunch().get_reference_particle().get_repetition();
  trajectory_length
    = get_bunch().get_reference_particle().get_trajectory_length();
  for (int idxtrk=0; idxtrk<num_tracks; ++idxtrk) {
    Track_status *dtsptr = &diag_track_status[idxtrk];
    if (dtsptr->found || first_search) {
      int index;
      dtsptr->found = false;
      if ((dtsptr->last_index > -1) &&
	  (dtsptr->last_index < get_bunch().get_local_num())) {
	if (dtsptr->particle_id
	    == static_cast<int > (get_bunch().get_local_particles()[Bunch::id][dtsptr->last_index])) {
	  index = dtsptr->last_index;
	  dtsptr->found = true;
	}
      }
      if (!(dtsptr->found)) {
	index = 0;
	while ((index < get_bunch().get_local_num())
	       && (dtsptr->particle_id
		   != static_cast<int > (get_bunch().get_local_particles()[index][Bunch::id]))) {
	  index += 1;
	}
	if (index < get_bunch().get_local_num()) {
	  dtsptr->found = true;
	  dtsptr->last_index = index;
	} else {
	  dtsptr->found = false;
	}
      }
      if (dtsptr->found) {
	if (first_search) {
	  get_write_helper();
	}
	track_coords[idxtrk][0] = get_bunch().get_local_particles()[index][0];
	track_coords[idxtrk][1] = get_bunch().get_local_particles()[index][1];
	track_coords[idxtrk][2] = get_bunch().get_local_particles()[index][2];
	track_coords[idxtrk][3] = get_bunch().get_local_particles()[index][3];
	track_coords[idxtrk][4] = get_bunch().get_local_particles()[index][4];
	track_coords[idxtrk][5] = get_bunch().get_local_particles()[index][5];
	track_coords[idxtrk][6] = get_bunch().get_local_particles()[index][6];
      } else {
	track_coords[idxtrk][0] = 0.0;
	track_coords[idxtrk][1] = 0.0;
	track_coords[idxtrk][2] = 0.0;
	track_coords[idxtrk][3] = 0.0;
	track_coords[idxtrk][4] = 0.0;
	track_coords[idxtrk][5] = 0.0;
	track_coords[idxtrk][6] = -static_cast <double>(dtsptr->particle_id);
      }	
      s = get_bunch().get_reference_particle().get_s();
      repetition = get_bunch().get_reference_particle().get_repetition();
      trajectory_length
	= get_bunch().get_reference_particle().get_trajectory_length();
      first_search = false;
    }
  }
}

void
Diagnostics_bulk_track::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_coords = new Hdf5_serial_writer<MArray2d_ref > (file_sptr,
                "track_coords");
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");
        writer_repetition = new Hdf5_serial_writer<int > (file_sptr,
                "repetition");
        writer_trajectory_length = new Hdf5_serial_writer<double > (file_sptr,
                "trajectory_length");
        have_writers = true;
    }
}

void
Diagnostics_bulk_track::write()
{
    const int max_writers = 4; // for now, may be dynamic in future
    int num_cycles = (get_bunch().get_comm().get_size() + max_writers - 1) / max_writers;
    Logger iocclog("iocycle_bulk_track", iocc_verbose);

    for (int cycle = 0; cycle<num_cycles; ++cycle) {
	iocclog << "start cycle " << cycle << std::endl;
        int cycle_min = cycle * max_writers;
        int cycle_max = (cycle + 1) * max_writers;
        if ((get_bunch().get_comm().get_rank() >= cycle_min)
	    && (get_bunch().get_comm().get_rank() < cycle_max)) {
            iocclog << "start write" << std::endl;

	    get_bunch().convert_to_state(get_bunch().fixed_z_lab);
	    init_writers(get_write_helper().get_hdf5_file_sptr());
	    writer_coords->append(track_coords);
	    writer_s->append(s);
	    writer_repetition->append(repetition);
	    writer_trajectory_length->append(trajectory_length);
	    get_write_helper().finish_write();
	}
	iocclog << "end write" << std::endl;
	MPI_Barrier(get_bunch().get_comm().get());
    }
}

template<class Archive>
    void
    Diagnostics_bulk_track::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
	    & BOOST_SERIALIZATION_NVP(num_tracks)
	    & BOOST_SERIALIZATION_NVP(iocc_verbose)
                & BOOST_SERIALIZATION_NVP(have_writers)
                & BOOST_SERIALIZATION_NVP(first_search)
	  & BOOST_SERIALIZATION_NVP(diag_track_status)
                & BOOST_SERIALIZATION_NVP(s)
                & BOOST_SERIALIZATION_NVP(writer_s)
                & BOOST_SERIALIZATION_NVP(repetition)
                & BOOST_SERIALIZATION_NVP(writer_repetition)
                & BOOST_SERIALIZATION_NVP(trajectory_length)
                & BOOST_SERIALIZATION_NVP(writer_trajectory_length)
                & BOOST_SERIALIZATION_NVP(track_coords)
                & BOOST_SERIALIZATION_NVP(writer_coords);
    }

template<class Archive>
    void
    Diagnostics_bulk_track::Track_status::serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP(found) &
	BOOST_SERIALIZATION_NVP(last_index) &
	BOOST_SERIALIZATION_NVP(particle_id);
    }

template
void
Diagnostics_bulk_track::Track_status::serialize<boost::archive::binary_oarchive >(
							       boost::archive::binary_oarchive &ar, const unsigned int version);

template
void
Diagnostics_bulk_track::Track_status::serialize<boost::archive::xml_oarchive >(
							       boost::archive::xml_oarchive &ar, const unsigned int version);

template
void
Diagnostics_bulk_track::Track_status::serialize<boost::archive::binary_iarchive >(
							       boost::archive::binary_iarchive &ar, const unsigned int version);

template
void
Diagnostics_bulk_track::Track_status::serialize<boost::archive::xml_iarchive >(
							       boost::archive::xml_iarchive &ar, const unsigned int version);


template
void
Diagnostics_bulk_track::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_bulk_track::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_bulk_track::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_bulk_track::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_bulk_track::~Diagnostics_bulk_track()
{
    if (have_writers) {
        delete writer_trajectory_length;
        delete writer_repetition;
        delete writer_s;
        delete writer_coords;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_bulk_track)

