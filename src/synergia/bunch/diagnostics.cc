#include "diagnostics.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "synergia/utils/eigen2/Eigen/Core"
#include "synergia/utils/eigen2/Eigen/LU"
#include <stdexcept>

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

MArray1d
Diagnostics::calculate_mean(Bunch const& bunch)
{
    MArray1d mean(boost::extents[6]);
    double sum[6] = { 0, 0, 0, 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            sum[i] += particles[part][i];
        }
    }
    MPI_Allreduce(sum, mean.origin(), 6, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 6; ++i) {
        mean[i] /= bunch.get_total_num();
    }
    return mean;
}

MArray1d
Diagnostics::calculate_std(Bunch const& bunch, MArray1d_ref const& mean)
{
    MArray1d std(boost::extents[6]);
    double sum[6] = { 0, 0, 0, 0, 0, 0 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            double diff = particles[part][i] - mean[i];
            sum[i] += diff * diff;
        }
    }
    MPI_Allreduce(sum, std.origin(), 6, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 6; ++i) {
        std[i] = std::sqrt(std[i] / bunch.get_total_num());
    }
    return std;
}

MArray2d
Diagnostics::calculate_mom2(Bunch const& bunch, MArray1d_ref const& mean)
{
    MArray2d mom2(boost::extents[6][6]);
    MArray2d sum2(boost::extents[6][6]);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            sum2[i][j] = 0.0;
        }
    }
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            double diff_i = particles[part][i] - mean[i];
            for (int j = 0; j <= i; ++j) {
                double diff_j = particles[part][j] - mean[j];
                sum2[i][j] += diff_i * diff_j;
            }
        }
    }
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            sum2[i][j] = sum2[j][i];
        }
    }
    MPI_Allreduce(sum2.origin(), mom2.origin(), 36, MPI_DOUBLE, MPI_SUM,
            bunch.get_comm().get());
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            mom2[i][j] = mom2[j][i] = mom2[i][j] / bunch.get_total_num();
        }
    }
    return mom2;
}

MArray1d
Diagnostics::calculate_bunchmin(Bunch const& bunch)
{
    MArray1d bunchmin(boost::extents[3]);
    double lmin[3] = { 1.0e100, 1.0e100, 1.0e100 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        if (particles[part][0] < lmin[0]) {
            lmin[0] = particles[part][0];
        }
        if (particles[part][2] < lmin[1]) {
            lmin[1] = particles[part][2];
        }
        if (particles[part][4] < lmin[2]) {
            lmin[2] = particles[part][4];
        }

    }
    MPI_Allreduce(lmin, bunchmin.origin(), 3, MPI_DOUBLE, MPI_MIN,
            bunch.get_comm().get());

    return bunchmin;
}

MArray1d
Diagnostics::calculate_bunchmax(Bunch const& bunch)
{
    MArray1d bunchmax(boost::extents[3]);
    double lmax[3] = { -1.0e100, -1.0e100, -1.0e100 };
    Const_MArray2d_ref particles(bunch.get_local_particles());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        if (particles[part][0] > lmax[0]) {
            lmax[0] = particles[part][0];
        }
        if (particles[part][2] > lmax[1]) {
            lmax[1] = particles[part][2];
        }
        if (particles[part][4] > lmax[2]) {
            lmax[2] = particles[part][4];
        }

    }
    MPI_Allreduce(lmax, bunchmax.origin(), 3, MPI_DOUBLE, MPI_MAX,
            bunch.get_comm().get());

    return bunchmax;
}

void
Diagnostics::update_and_write()
{
    update();
    write();
}

Diagnostics_basic::Diagnostics_basic(Bunch_sptr bunch_sptr,
        std::string const& filename) :
    bunch_sptr(bunch_sptr), filename(filename), have_writers(false),
            mean(boost::extents[6]), std(boost::extents[6]),
            write_helper(filename, true, bunch_sptr->get_comm())
{
}

bool
Diagnostics_basic::is_serial() const
{
    return true;
}

void
Diagnostics_basic::update()
{
    s = bunch_sptr->get_reference_particle().get_s();
    repetition = bunch_sptr->get_reference_particle().get_repetition();
    trajectory_length
            = bunch_sptr->get_reference_particle().get_trajectory_length();
    num_particles = bunch_sptr->get_total_num();
    real_num_particles = bunch_sptr->get_real_num();
    mean = Diagnostics::calculate_mean(*bunch_sptr);
    std = Diagnostics::calculate_std(*bunch_sptr, mean);
}

double
Diagnostics_basic::get_s() const
{
    return s;
}

int
Diagnostics_basic::get_repetition() const
{
    return repetition;
}

double
Diagnostics_basic::get_trajectory_length() const
{
    return trajectory_length;
}

int
Diagnostics_basic::get_num_particles() const
{
    return num_particles;
}

double
Diagnostics_basic::get_real_num_particles() const
{
    return real_num_particles;
}

Const_MArray1d_ref
Diagnostics_basic::get_mean() const
{
    return mean;
}

Const_MArray1d_ref
Diagnostics_basic::get_std() const
{
    return std;
}

const MArray1d
Diagnostics_basic::get_bunchmin() const
{
    MArray1d bunchmin;
    bunchmin = Diagnostics::calculate_bunchmin(*bunch_sptr);
    return bunchmin;
}

const MArray1d
Diagnostics_basic::get_bunchmax() const
{
    MArray1d bunchmax;
    bunchmax = Diagnostics::calculate_bunchmax(*bunch_sptr);
    return bunchmax;
}

void
Diagnostics_basic::init_writers(H5::H5File & file)
{
    if (!have_writers) {
        writer_s = new Hdf5_serial_writer<double > (file, "s");
        writer_repetition = new Hdf5_serial_writer<int > (file, "repetition");
        writer_trajectory_length = new Hdf5_serial_writer<double > (file,
                "trajectory_length");
        writer_num_particles = new Hdf5_serial_writer<int > (file,
                "num_particles");
        writer_real_num_particles = new Hdf5_serial_writer<double > (file,
                "real_num_particles");
        writer_mean = new Hdf5_serial_writer<MArray1d_ref > (file, "mean");
        writer_std = new Hdf5_serial_writer<MArray1d_ref > (file, "std");
        have_writers = true;
    }
}

void
Diagnostics_basic::write()
{
    if (write_helper.write_locally()) {
        init_writers(write_helper.get_file());
        writer_s->append(s);
        writer_repetition->append(repetition);
        writer_trajectory_length->append(trajectory_length);
        writer_num_particles->append(num_particles);
        writer_real_num_particles->append(real_num_particles);
        writer_mean->append(mean);
        writer_std->append(std);
        write_helper.finish_write();
    }
}

Diagnostics_basic::~Diagnostics_basic()
{
    if (have_writers) {
        delete writer_std;
        delete writer_mean;
        delete writer_real_num_particles;
        delete writer_num_particles;
        delete writer_trajectory_length;
        delete writer_repetition;
        delete writer_s;
    }
}

void
Diagnostics_full2::update_full2()
{
    MArray2d sum2(boost::extents[6][6]);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            sum2[i][j] = 0.0;
        }
    }
    Const_MArray2d_ref particles(bunch_sptr->get_local_particles());
    for (int part = 0; part < bunch_sptr->get_local_num(); ++part) {
        for (int i = 0; i < 6; ++i) {
            double diff_i = particles[part][i] - mean[i];
            for (int j = 0; j <= i; ++j) {
                double diff_j = particles[part][j] - mean[j];
                sum2[i][j] += diff_i * diff_j;
            }
        }
    }
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            sum2[i][j] = sum2[j][i];
        }
    }
    MPI_Allreduce(sum2.origin(), mom2.origin(), 36, MPI_DOUBLE, MPI_SUM,
            bunch_sptr->get_comm().get());
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            mom2[i][j] = mom2[j][i] = mom2[i][j] / bunch_sptr->get_total_num();
        }
        std[i] = std::sqrt(mom2[i][i]);
    }
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            corr[i][j] = corr[j][i] = mom2[i][j] / std::sqrt(
                    mom2[i][i] * mom2[j][j]);
        }
    }
}

void
Diagnostics_full2::update_emittances()
{
    Matrix<double, 6, 6 > mom2_matrix(mom2.origin());
    emitx = std::sqrt(
            mom2_matrix.block<2, 2 > (Bunch::x, Bunch::x).determinant());
    emity = std::sqrt(
            mom2_matrix.block<2, 2 > (Bunch::y, Bunch::y).determinant());
    emitz = std::sqrt(
            mom2_matrix.block<2, 2 > (Bunch::z, Bunch::z).determinant());
    emitxy = std::sqrt(
            mom2_matrix.block<4, 4 > (Bunch::x, Bunch::x).determinant());
    emitxyz = std::sqrt(mom2_matrix.determinant());
}

Diagnostics_full2::Diagnostics_full2(Bunch_sptr bunch_sptr,
        std::string const& filename) :
    bunch_sptr(bunch_sptr), filename(filename), have_writers(false),
            mean(boost::extents[6]), std(boost::extents[6]),
            mom2(boost::extents[6][6]), corr(boost::extents[6][6]),
            write_helper(filename, true, bunch_sptr->get_comm())
{          
}

bool
Diagnostics_full2::is_serial() const
{
    return true;
}

void
Diagnostics_full2::update()
{   
    bunch_sptr->convert_to_state(bunch_sptr->fixed_z);
    s = bunch_sptr->get_reference_particle().get_s();
    repetition = bunch_sptr->get_reference_particle().get_repetition();
    trajectory_length
            = bunch_sptr->get_reference_particle().get_trajectory_length();
    num_particles = bunch_sptr->get_total_num();
    real_num_particles = bunch_sptr->get_real_num();
    mean = calculate_mean(*bunch_sptr);
    update_full2();
    update_emittances();
}

double
Diagnostics_full2::get_s() const
{
    return s;
}

int
Diagnostics_full2::get_repetition() const
{
    return repetition;
}

double
Diagnostics_full2::get_trajectory_length() const
{
    return trajectory_length;
}

int
Diagnostics_full2::get_num_particles() const
{
    return num_particles;
}

double
Diagnostics_full2::get_real_num_particles() const
{
    return real_num_particles;
}

Const_MArray1d_ref
Diagnostics_full2::get_mean() const
{
    return mean;
}

Const_MArray1d_ref
Diagnostics_full2::get_std() const
{
    return std;
}

Const_MArray2d_ref
Diagnostics_full2::get_mom2() const
{
    return mom2;
}

Const_MArray2d_ref
Diagnostics_full2::get_corr() const
{
    return corr;
}

double
Diagnostics_full2::get_emitx() const
{
    return emitx;
}

double
Diagnostics_full2::get_emity() const
{
    return emity;
}

double
Diagnostics_full2::get_emitz() const
{
    return emitz;
}

double
Diagnostics_full2::get_emitxy() const
{
    return emitxy;
}

double
Diagnostics_full2::get_emitxyz() const
{
    return emitxyz;
}

void
Diagnostics_full2::init_writers(H5::H5File & file)
{
    if (!have_writers) {
        writer_s = new Hdf5_serial_writer<double > (file, "s");
        writer_repetition = new Hdf5_serial_writer<int > (file, "repetition");
        writer_trajectory_length = new Hdf5_serial_writer<double > (file,
                "trajectory_length");
        writer_num_particles = new Hdf5_serial_writer<int > (file,
                "num_particles");
        writer_real_num_particles = new Hdf5_serial_writer<double > (file,
                "real_num_particles");
        writer_mean = new Hdf5_serial_writer<MArray1d_ref > (file, "mean");
        writer_std = new Hdf5_serial_writer<MArray1d_ref > (file, "std");
        writer_mom2 = new Hdf5_serial_writer<MArray2d_ref > (file, "mom2");
        writer_corr = new Hdf5_serial_writer<MArray2d_ref > (file, "corr");
        writer_emitx = new Hdf5_serial_writer<double > (file, "emitx");
        writer_emity = new Hdf5_serial_writer<double > (file, "emity");
        writer_emitz = new Hdf5_serial_writer<double > (file, "emitz");
        writer_emitxy = new Hdf5_serial_writer<double > (file, "emitxy");
        writer_emitxyz = new Hdf5_serial_writer<double > (file, "emitxyz");
        have_writers = true;
    }
}

void
Diagnostics_full2::write()
{
    if (write_helper.write_locally()) {
        init_writers(write_helper.get_file());
        writer_s->append(s);
        writer_repetition->append(repetition);
        writer_trajectory_length->append(trajectory_length);
        writer_num_particles->append(num_particles);
        writer_real_num_particles->append(real_num_particles);
        writer_mean->append(mean);
        writer_std->append(std);
        writer_mom2->append(mom2);
        writer_corr->append(corr);
        writer_emitx->append(emitx);
        writer_emity->append(emity);
        writer_emitz->append(emitz);
        writer_emitxy->append(emitxy);
        writer_emitxyz->append(emitxyz);
        write_helper.finish_write();
    }
}

Diagnostics_full2::~Diagnostics_full2()
{
    if (have_writers) {
        delete writer_emitxyz;
        delete writer_emitxy;
        delete writer_emitz;
        delete writer_emity;
        delete writer_emitx;
        delete writer_corr;
        delete writer_mom2;
        delete writer_std;
        delete writer_mean;
        delete writer_real_num_particles;
        delete writer_num_particles;
        delete writer_trajectory_length;
        delete writer_repetition;
        delete writer_s;
    }
}

Diagnostics_particles::Diagnostics_particles(Bunch_sptr bunch_sptr,
        std::string const& filename,  int max_particles) :
    bunch_sptr(bunch_sptr), filename(filename), max_particles(max_particles),
            have_writers(false), write_helper(filename, false,
                    bunch_sptr->get_comm())
{
}

Diagnostics_particles::Diagnostics_particles(Bunch_sptr bunch_sptr,
        std::string const& filename, int write_skip, int max_particles) :
    bunch_sptr(bunch_sptr), filename(filename), max_particles(max_particles),
            have_writers(false), write_helper(filename, false, write_skip,
                    bunch_sptr->get_comm())
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
            int particle_id = particles[part][Bunch::id];
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
        std::vector<int > const& local_nums, H5::H5File & file)
{
    int myrank = bunch_sptr->get_comm().get_rank();
    int size = bunch_sptr->get_comm().get_size();
    Hdf5_chunked_array2d_writer
            writer_particles(
                    file,
                    "particles",
                    bunch_sptr->get_local_particles()[boost::indices[range(0, 1)][range()]]);
    for (int rank = 0; rank < size; ++rank) {
        int local_num = local_nums[rank];
        if (rank == myrank) {
            write_selected_particles(writer_particles,
                    bunch_sptr->get_local_particles(), local_num,
                    min_particle_id, max_particle_id);
        } else {
            MPI_Status status;
            MArray2d received(boost::extents[local_num][7]);
            int message_size = 7 * local_num;
            MPI_Comm comm = bunch_sptr->get_comm().get();
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
    int local_num = bunch_sptr->get_local_num();
    void * send_buffer =
            (void*) bunch_sptr->get_local_particles()[boost::indices[range(0,
                    local_num)][range()]].origin();
    int status;
    int message_size = 7 * local_num;
    int receiver = write_helper.get_writer_rank();
    int rank = bunch_sptr->get_comm().get_rank();
    MPI_Comm comm = bunch_sptr->get_comm().get();
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

    int writer_rank= write_helper.get_writer_rank();
    MPI_Comm comm = bunch_sptr->get_comm().get();
    int icount;
    icount=write_helper.get_count();
    MPI_Bcast ((void *) &icount, 1, MPI_INT, writer_rank, comm );


//    std::cout<<" icount ="<<icount<<"  count="<< write_helper.get_count() <<" on rank ="<< bunch_sptr->get_comm().get_rank()<<std::endl;
     
     
    if (icount % write_helper.get_iwrite_skip() !=0 ) 
    {   
         if (write_helper.write_locally()) write_helper.increment_count();
         return;
    }
   
  
      
   
    int local_num = bunch_sptr->get_local_num();
    int num_procs = bunch_sptr->get_comm().get_size();
    std::vector<int > local_nums(num_procs);
    void * local_nums_buf = (void *) &local_nums[0];
    int root = write_helper.get_writer_rank();
    int status;
    status = MPI_Gather((void*) &local_num, 1, MPI_INT, local_nums_buf, 1,
            MPI_INT, root, comm);
    if (status != MPI_SUCCESS) {
        throw std::runtime_error(
                "Diagnostics_particles::write: MPI_Gather failed.");
    }

    if (write_helper.write_locally()) {
        H5::H5File file = write_helper.get_file();
        receive_other_local_particles(local_nums, file);

        Hdf5_writer<double > writer_pz(file, "pz");
        double pz = bunch_sptr->get_reference_particle().get_momentum();
        writer_pz.write(pz);
        Hdf5_writer<double > writer_tlen(file, "tlen");
        double tlen =
                bunch_sptr->get_reference_particle().get_trajectory_length();
        writer_tlen.write(tlen);
        Hdf5_writer<int > writer_rep(file, "rep");
        int rep = bunch_sptr->get_reference_particle().get_repetition();
        writer_rep.write(rep);
        Hdf5_writer<double > writer_s(file, "s");
        double s = bunch_sptr->get_reference_particle().get_s();
        writer_s.write(s);
        int local_num = bunch_sptr->get_local_num();
        write_helper.finish_write();
    } else {
        send_local_particles();
    }
}

Diagnostics_particles::~Diagnostics_particles()
{
}

Diagnostics_track::Diagnostics_track(Bunch_sptr bunch_sptr,
        std::string const& filename, int particle_id) :
    bunch_sptr(bunch_sptr), filename(filename), have_writers(false),
            coords(boost::extents[6]), found(false), first_search(true),
            particle_id(particle_id), last_index(-1)
{
}

bool
Diagnostics_track::is_serial() const
{
    return true;
}

void
Diagnostics_track::update()
{
    repetition = bunch_sptr->get_reference_particle().get_repetition();
    trajectory_length
            = bunch_sptr->get_reference_particle().get_trajectory_length();
    if (found || first_search) {
        int index;
        found = false;
        if ((last_index > -1) && (last_index < bunch_sptr->get_local_num())) {
            if (particle_id
                    == static_cast<int > (bunch_sptr->get_local_particles()[Bunch::id][last_index])) {
                index = last_index;
                found = true;
            }
        }
        if (!found) {
            index = 0;
            while ((index < bunch_sptr->get_local_num())
                    && (particle_id
                            != static_cast<int > (bunch_sptr->get_local_particles()[index][Bunch::id]))) {
                index += 1;
            }
            if (index < bunch_sptr->get_local_num()) {
                found = true;
            } else {
                found = false;
            }
        }
        if (found) {
            if (first_search) {
                write_helper_sptr = Diagnostics_write_helper_sptr(
                        new Diagnostics_write_helper(filename, true,
                                bunch_sptr->get_comm(),
                                bunch_sptr->get_comm().get_rank()));
            }
            coords[0] = bunch_sptr->get_local_particles()[index][0];
            coords[1] = bunch_sptr->get_local_particles()[index][1];
            coords[2] = bunch_sptr->get_local_particles()[index][2];
            coords[3] = bunch_sptr->get_local_particles()[index][3];
            coords[4] = bunch_sptr->get_local_particles()[index][4];
            coords[5] = bunch_sptr->get_local_particles()[index][5];
            s = bunch_sptr->get_reference_particle().get_s();
            repetition = bunch_sptr->get_reference_particle().get_repetition();
            trajectory_length
                    = bunch_sptr->get_reference_particle().get_trajectory_length();
        }
        first_search = false;
    }
}

void
Diagnostics_track::init_writers(H5::H5File & file)
{
    if (!have_writers) {
        writer_coords = new Hdf5_serial_writer<MArray1d_ref > (file, "coords");
        writer_s = new Hdf5_serial_writer<double > (file, "s");
        writer_repetition = new Hdf5_serial_writer<int > (file, "repetition");
        writer_trajectory_length = new Hdf5_serial_writer<double > (file,
                "trajectory_length");
        have_writers = true;
    }
}

void
Diagnostics_track::write()
{

    if (found) {
        init_writers(write_helper_sptr->get_file());
        writer_coords->append(coords);
        writer_s->append(s);
        writer_repetition->append(repetition);
        writer_trajectory_length->append(trajectory_length);
        write_helper_sptr->finish_write();
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
