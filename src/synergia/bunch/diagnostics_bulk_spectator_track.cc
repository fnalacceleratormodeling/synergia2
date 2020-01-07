#include <string>
#include <iostream>
#include "diagnostics_bulk_spectator_track.h"
#include "synergia/utils/parallel_utils.h"

const char Diagnostics_bulk_spectator_track::name[] = "diagnostics_bulk_spectator_track";

const int no_local_num_tracks = -1;

Diagnostics_bulk_spectator_track::Diagnostics_bulk_spectator_track(
        std::string const& filename,
        int num_tracks, 
        int offset, 
        std::string const& local_dir) 
    : Diagnostics(Diagnostics_bulk_spectator_track::name, filename, local_dir)
    , num_tracks(num_tracks)
    , local_num_tracks(no_local_num_tracks)
    , offset(offset)
    , local_offset(no_local_num_tracks)
    , have_writers(false)
    , first_search(true)
    , diag_track_status()
    , writer_s_n(0)
    , writer_repetition(0)
    , writer_s(0)
    , writer_pz(0)
    , track_coords(boost::extents[1][7])
    , writer_coords(0)
{
}

Diagnostics_bulk_spectator_track::Diagnostics_bulk_spectator_track() : have_writers(false)
{
}

bool
Diagnostics_bulk_spectator_track::is_serial() const
{
    return true;
}

//Diagnostics_write_helper *
//Diagnostics_bulk_spectator_track::new_write_helper_ptr()
//{
//    delete_write_helper_ptr();
//    std::stringstream sstream;
//    std::string filename(get_filename());
//    int idx = filename.rfind('.');
//    std::string filename_base, filename_suffix;
//    if (idx == std::string::npos) {
//        filename_base = filename;
//        filename_suffix = "";
//    } else {
//        filename_base = filename.substr(0, idx);
//        filename_suffix = filename.substr(idx);
//    }
//
//    sstream << filename_base;
//    sstream << "_";
//    sstream << std::setw(4);
//
//    sstream << std::setfill('0');
//    sstream << get_bunch().get_comm().get_rank();
//    sstream << filename_suffix;
//    return new Diagnostics_write_helper(sstream.str(), true,
//            get_bunch().get_comm(), get_bunch().get_comm().get_rank());
//}

void
Diagnostics_bulk_spectator_track::update()
{
    if (get_bunch().get_comm().has_this_rank())
    {
        if (local_num_tracks == no_local_num_tracks) 
        {
            local_num_tracks = decompose_1d_local(get_bunch().get_comm(), num_tracks);
            local_offset = decompose_1d_local(get_bunch().get_comm(), offset);
            track_coords.resize(boost::extents[local_num_tracks][7]);
        }

        // first time, set the status of the tracks that are being followed
        if (diag_track_status.empty()) 
        {
            if (local_num_tracks + local_offset > get_bunch().get_local_spectator_num()) 
            {
                local_num_tracks = get_bunch().get_local_spectator_num() - local_offset;
            }

            // loop over my local tracks
            for (int idxtrk = 0; idxtrk < local_num_tracks; ++idxtrk) 
            {
               // here's one
                Track_status dts;
                dts.found = true;
                dts.last_index = idxtrk+local_offset;
                dts.particle_id =
                        static_cast<int >(get_bunch().get_local_spectator_particles()[idxtrk+local_offset][Bunch::id]);
                diag_track_status.push_back(dts);
            }
        }

        get_bunch().convert_to_state(get_bunch().fixed_z_lab);
        repetition = get_bunch().get_reference_particle().get_repetition();
        s = get_bunch().get_reference_particle().get_s();

        // for the tracks I'm following, see if I know where they are
        for (int idxtrk = 0; idxtrk < local_num_tracks; ++idxtrk) 
        {
            // diag_track_status[0..local_num_tracks] has status for tracks originally
            // at local_particles[local_offset..local_offset+local_num_tracks]
            Track_status *dtsptr = &diag_track_status[idxtrk];
            if (dtsptr->found || first_search) 
            {
                int index = 0;
                dtsptr->found = false;

                if ((dtsptr->last_index > -1)
                        && (dtsptr->last_index < get_bunch().get_local_spectator_num())) 
                {
                    if (dtsptr->particle_id
                            == static_cast<int >(get_bunch().get_local_spectator_particles()[dtsptr->last_index][Bunch::id])) 
                    {
                        // this track is still where it was last time I looked
                        index = dtsptr->last_index;
                        dtsptr->found = true;
                    }
                }

                if (!(dtsptr->found)) 
                {
                    // it wasn't there, I have to look for it again
                    index = 0;
                    while ((index < get_bunch().get_local_spectator_num())
                           && (dtsptr->particle_id
                               != static_cast<int >(get_bunch().get_local_spectator_particles()[index][Bunch::id]))) 
                    {
                        index += 1;
                    }

                    if (index < get_bunch().get_local_spectator_num()) 
                    {
                        // found it
                        dtsptr->found = true;
                        dtsptr->last_index = index;
                    } 
                    else 
                    {
                        // It's not there'
                        dtsptr->found = false;
                    }
                }

                if (dtsptr->found) 
                {
                    if (first_search) 
                    {
                        get_write_helper();
                    }

                    // track_coords is parallel  to diag_track_status
                    // diag_track_status[0..local_num_tracks] has status for tracks originally
                    // at local_particles[local_offset..local_offset+local_num_tracks]
                    track_coords[ boost::indices[idxtrk][range()] ] =
                            get_bunch().get_local_spectator_particles()[boost::indices[index][range()]];
               } 
               else 
               {
                    track_coords[idxtrk][0] = 0.0;
                    track_coords[idxtrk][1] = 0.0;
                    track_coords[idxtrk][2] = 0.0;
                    track_coords[idxtrk][3] = 0.0;
                    track_coords[idxtrk][4] = 0.0;
                    track_coords[idxtrk][5] = 0.0;
                    track_coords[idxtrk][6] = -static_cast<double >(dtsptr->particle_id);
                }

                s_n = get_bunch().get_reference_particle().get_s_n();
                repetition = get_bunch().get_reference_particle().get_repetition();
                s = get_bunch().get_reference_particle().get_s();
                pz = get_bunch().get_reference_particle().get_momentum();
                first_search = false;
            }
        }
    }
}

void
Diagnostics_bulk_spectator_track::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );
        int chg = get_bunch().get_reference_particle().get_charge();
        file_sptr->write(chg, "charge");
        double pmass = fourp.get_mass();
        file_sptr->write(pmass, "mass");

        writer_coords = new Hdf5_serial_writer<MArray2d_ref >(file_sptr,
                                                              "track_coords");
        writer_s_n = new Hdf5_serial_writer<double >(file_sptr, "s_n");
        writer_repetition = new Hdf5_serial_writer<int >(file_sptr,
                                                         "repetition");
        writer_s = new Hdf5_serial_writer<double >(file_sptr,
                                                   "s");
        writer_pz = new Hdf5_serial_writer<double >(file_sptr, "pz");
        have_writers = true;
    }
}

void
Diagnostics_bulk_spectator_track::receive_other_local_coords(
        std::vector<int > const& local_nums)
{
    if (get_bunch().get_comm().has_this_rank()){
        int myrank = get_bunch().get_comm().get_rank();
        int size = get_bunch().get_comm().get_size();
        MArray2d all_coords(boost::extents[num_tracks][7]);
        int array_offset = 0;
        for (int rank = 0; rank < size; ++rank) {
            int local_num = local_nums[rank];
            if (rank == myrank) {
                for (int i = array_offset; i < array_offset + local_num; ++i) {
                    for (int j = 0; j < 7; ++j) {
                        all_coords[i][j] = track_coords[i - array_offset][j];
                    }
                }
            } else {
                MPI_Status status;
                int message_size = 7 * local_num;
                MPI_Comm comm = get_bunch().get_comm().get();
                int error = MPI_Recv((void*) &all_coords[array_offset][0],
                                     message_size, MPI_DOUBLE, rank, rank, comm, &status);
                if (error != MPI_SUCCESS) {
                    throw std::runtime_error(
                                "Diagnostics_bulk_spectator_track::receive_other_local_coords: MPI_Recv failed.");
                }
                array_offset += local_num;
            }
        }
        writer_coords->append(all_coords);
    }
}

void
Diagnostics_bulk_spectator_track::send_local_coords()
{
    if (get_bunch().get_comm().has_this_rank()){
        void * send_buffer =
                (void*) track_coords.origin();
        int status;
        int message_size = 7 * local_num_tracks;
        int receiver = get_write_helper().get_writer_rank();
        int rank = get_bunch().get_comm().get_rank();
        MPI_Comm comm = get_bunch().get_comm().get();
        status = MPI_Send(send_buffer, message_size, MPI_DOUBLE, receiver, rank,
                          comm);
        if (status != MPI_SUCCESS) {
            throw std::runtime_error(
                        "Diagnostics_bulk_spectator_track::send_local_coords: MPI_Send failed.");
        }
    }
}

void
Diagnostics_bulk_spectator_track::write()
{   if (get_bunch().get_comm().has_this_rank()){
        MPI_Comm comm = get_bunch().get_comm().get();
        int num_procs = get_bunch().get_comm().get_size();

        std::vector<int > all_local_num_tracks(num_procs);
        void * local_num_tracks_buf = (void *) &all_local_num_tracks[0];
        int root = get_write_helper().get_writer_rank();
        int status;
        status = MPI_Gather((void*) &local_num_tracks, 1, MPI_INT,
                            local_num_tracks_buf, 1, MPI_INT, root, comm);
        if (status != MPI_SUCCESS) {
            throw std::runtime_error(
                        "Diagnostics_bulk_spectator_track::write: MPI_Gather failed.");
        }

        if (get_write_helper().write_locally()) {
            init_writers(get_write_helper().get_hdf5_file_sptr());
            receive_other_local_coords(all_local_num_tracks);
            writer_s_n->append(s_n);
            writer_repetition->append(repetition);
            writer_s->append(s);
            writer_pz->append(pz);
            get_write_helper().finish_write();
        } else {
            send_local_coords();
        }
    }
}

template<class Archive>
void
Diagnostics_bulk_spectator_track::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
    ar & BOOST_SERIALIZATION_NVP(num_tracks);
    ar & BOOST_SERIALIZATION_NVP(local_num_tracks);
    ar & BOOST_SERIALIZATION_NVP(offset);
    ar & BOOST_SERIALIZATION_NVP(local_offset);
    ar & BOOST_SERIALIZATION_NVP(have_writers);
    ar & BOOST_SERIALIZATION_NVP(first_search);
    ar & BOOST_SERIALIZATION_NVP(diag_track_status);
    ar & BOOST_SERIALIZATION_NVP(s_n);
    ar & BOOST_SERIALIZATION_NVP(writer_s_n);
    ar & BOOST_SERIALIZATION_NVP(repetition);
    ar & BOOST_SERIALIZATION_NVP(writer_repetition);
    ar & BOOST_SERIALIZATION_NVP(s);
    ar & BOOST_SERIALIZATION_NVP(writer_s);
    ar & BOOST_SERIALIZATION_NVP(pz);
    ar & BOOST_SERIALIZATION_NVP(writer_pz);
    ar & BOOST_SERIALIZATION_NVP(track_coords);
    ar & BOOST_SERIALIZATION_NVP(writer_coords);
}

template<class Archive>
void
Diagnostics_bulk_spectator_track::Track_status::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(found) &
            BOOST_SERIALIZATION_NVP(last_index) &
            BOOST_SERIALIZATION_NVP(particle_id);
}

template
void
Diagnostics_bulk_spectator_track::Track_status::serialize<boost::archive::binary_oarchive >(
boost::archive::binary_oarchive &ar, const unsigned int version);

template
void
Diagnostics_bulk_spectator_track::Track_status::serialize<boost::archive::xml_oarchive >(
boost::archive::xml_oarchive &ar, const unsigned int version);

template
void
Diagnostics_bulk_spectator_track::Track_status::serialize<boost::archive::binary_iarchive >(
boost::archive::binary_iarchive &ar, const unsigned int version);

template
void
Diagnostics_bulk_spectator_track::Track_status::serialize<boost::archive::xml_iarchive >(
boost::archive::xml_iarchive &ar, const unsigned int version);


template
void
Diagnostics_bulk_spectator_track::serialize<boost::archive::binary_oarchive >(
boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_bulk_spectator_track::serialize<boost::archive::xml_oarchive >(
boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_bulk_spectator_track::serialize<boost::archive::binary_iarchive >(
boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_bulk_spectator_track::serialize<boost::archive::xml_iarchive >(
boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_bulk_spectator_track::~Diagnostics_bulk_spectator_track()
{
    if (have_writers) {
        delete writer_pz;
        delete writer_s;
        delete writer_repetition;
        delete writer_s_n;
        delete writer_coords;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_bulk_spectator_track)

