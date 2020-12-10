
#include <iomanip>

#include "synergia/bunch/bunch_particles.h"
#include "synergia/utils/parallel_utils.h"
#include "synergia/utils/hdf5_file.h"

namespace
{
    struct Particle_id_offset
    {
        static int offset;

        static int get(int request_num, Commxx const& comm)
        {
            MPI_Bcast((void *)&offset, 1, MPI_INT, 0, comm);
            int old_offset = offset;
            int total_num;
            MPI_Reduce((void*)&request_num, (void*)&total_num, 1, 
                    MPI_INT, MPI_SUM, 0, comm);
            offset += total_num;
            return old_offset;
        }
    };

    int Particle_id_offset::offset = 0;

    struct particle_id_assigner
    {
        Particles parts;
        int64_t offset;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { parts(i, 6) = i + offset; }
    };

    struct particle_masks_initializer
    {
        ParticleMasks masks;
        const int num;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { masks(i) = i<num ? 1 : 0; }
    };

    // Kokkos functors
    struct particle_copier_many
    {
        ConstParticles src;
        karray2d_row_dev dst;

        ConstParticleMasks s_masks;
        ParticleMasks      d_masks;

        int idx;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            dst(i, 0) = src(idx+i, 0);
            dst(i, 1) = src(idx+i, 1);
            dst(i, 2) = src(idx+i, 2);
            dst(i, 3) = src(idx+i, 3);
            dst(i, 4) = src(idx+i, 4);
            dst(i, 5) = src(idx+i, 5);
            dst(i, 6) = src(idx+i, 6);

            d_masks(i) = s_masks(idx+i);
        }
    };

    struct particle_copier_one
    {
        ConstParticles src;
        karray1d_row_dev dst;

        ConstParticleMasks s_masks;
        ParticleMasks d_masks;

        int idx;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            dst(0) = src(idx+i, 0);
            dst(1) = src(idx+i, 1);
            dst(2) = src(idx+i, 2);
            dst(3) = src(idx+i, 3);
            dst(4) = src(idx+i, 4);
            dst(5) = src(idx+i, 5);
            dst(6) = src(idx+i, 6);

            d_masks(0) = s_masks(idx+i);
        }
    };

    struct particle_id_checker
    {
        ConstParticles parts;
        int idx;
        int pid;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& match) const
        { match = (((int)parts(idx, 6)) == pid); }
    };

    struct particle_finder
    {
        ConstParticles parts;
        int pid;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& idx) const
        { if (((int)parts(i, 6)) == pid) idx = i; }
    };

    struct particle_zeroer
    {
        Particles parts;
        int offset;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { for (int j=0; j<7; ++j) parts(offset+i, j) = 0.0; }
    };

    struct mask_reducer
    {
        ConstParticleMasks masks;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& valid) const
        { if(masks(i)) ++valid; }
    };

    struct particle_injector
    {
        Particles dst;
        ParticleMasks dst_masks;

        ConstParticles src;
        ConstParticleMasks src_masks;

        karray1d_dev ref_st_diff;
        karray1d_dev tgt_st;
        karray1d_dev inj_st;

        int off;
        double pdiff;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            // space-like coordinates
            dst(off+i, 0) = src(i, 0);// + ref_st_diff(0);
            dst(off+i, 2) = src(i, 2);// + ref_st_diff(2);
            dst(off+i, 4) = src(i, 4);// + ref_st_diff(4);

            // npx and npy coordinates are scaled with p_ref which can be different
            // for different bunches
            dst(off+i, 1) = pdiff * (src(i, 1) - inj_st(1)) + tgt_st(1);
            dst(off+i, 3) = pdiff * (src(i, 3) - inj_st(3)) + tgt_st(3);

            // ndp coordinate is delta-p scaled with pref
            dst(off+i, 5) = pdiff * (1.0 + src(i, 5) - inj_st(5)) + tgt_st(5) - 1.0;

            // particle id
            dst(off+i, 6) = src(i, 6);

            dst_masks(off+i) = src_masks(i);
        }
    };

    struct fixed_z_to_t_converter
    {
        Particles parts;
        ConstParticleMasks masks;

        double p_ref;
        double beta;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            if (masks(i))
            {
                double p = p_ref + parts(i, 5) * p_ref;
                double px = parts(i, 1) * p_ref;
                double py = parts(i, 3) * p_ref;
                double pz2 = p*p - px*px - py*py;
                double pz = sqrt(pz2);

                parts(i, 4) = - parts(i, 4) * beta;
                parts(i, 5) = pz / p_ref;
            }
        }
    };

    struct fixed_t_to_z_converter
    {
        Particles parts;
        ConstParticleMasks masks;

        double p_ref;
        double beta;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            if (masks(i))
            {
                double px = parts(i, 1) * p_ref;
                double py = parts(i, 3) * p_ref;
                double pz = parts(i, 5) * p_ref;
                double p = sqrt(px*px + py*py + pz*pz);

                parts(i, 4) = - parts(i, 4) / beta;
                parts(i, 5) = (p-p_ref) / p_ref;
            }
        }
    };

    struct discarded_particle_mover
    {
        Kokkos::View<int*> counter;

        ConstParticles   parts;
        ParticleMasks    discards;
        karray2d_row_dev discarded_parts;

        discarded_particle_mover(
                ConstParticles   const& parts,
                ParticleMasks    const& discards,
                karray2d_row_dev const& discarded_parts )
            : counter("counter", 1)
            , parts(parts)
            , discards(discards)
            , discarded_parts(discarded_parts)
        { }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            if (discards(i))
            {
                int pos = Kokkos::atomic_fetch_add(counter.data(), 1);

                discarded_parts(pos, 0) = parts(i, 0);
                discarded_parts(pos, 1) = parts(i, 1);
                discarded_parts(pos, 2) = parts(i, 2);
                discarded_parts(pos, 3) = parts(i, 3);
                discarded_parts(pos, 4) = parts(i, 4);
                discarded_parts(pos, 5) = parts(i, 5);
                discarded_parts(pos, 6) = parts(i, 6);
            }
        }
    };
}


template<>
void bunch_particles_t<double>::assign_ids(int local_offset, Commxx const& comm)
{
    int request_num = (comm.rank() == 0) ? n_total : 0;
    int global_offset = Particle_id_offset::get(request_num, comm);

    particle_id_assigner pia{parts, local_offset + global_offset};
    Kokkos::parallel_for(n_active, pia);
}

template<>
bunch_particles_t<double>::bunch_particles_t( 
        ParticleGroup pg,
        int total, int reserved, 
        Commxx const& comm)
    : group(pg)
    , label(pg==PG::regular ? "particles" : "spectators")
    , n_valid(0)
    , n_active(0)
    , n_reserved(0)
    , n_total(total)
    , n_last_discarded(0)
    , poffset(0)
    , parts()
    , masks()
    , discards()
    , hparts(Kokkos::create_mirror_view(parts))
    , hmasks(Kokkos::create_mirror_view(masks))
    , hdiscards(Kokkos::create_mirror_view(discards))
{
    // minimum reserved
    if (reserved < total) reserved = total;

    if (!comm.is_null() && reserved) 
    {
        int mpi_size = comm.size();
        int mpi_rank = comm.rank();

        std::vector<int> offsets_t(mpi_size);
        std::vector<int> counts_t(mpi_size);
        decompose_1d(comm, total, offsets_t, counts_t);

        std::vector<int> offsets_r(mpi_size);
        std::vector<int> counts_r(mpi_size);
        decompose_1d(comm, reserved, offsets_r, counts_r);

        // local_num
        n_valid    = counts_t[mpi_rank];
        n_active   = counts_t[mpi_rank];
        n_reserved = counts_r[mpi_rank];

        if (n_active % gsv_t::size())
        {
            int padded = n_active + gsv_t::size() - n_active%gsv_t::size();
            if (n_reserved < padded) n_reserved = padded;
        }

        // local_num offset
        poffset = 0;
        for(int i=0; i<mpi_rank; ++i) 
            poffset+=counts_t[i];

        // allocate
#ifdef NO_PADDING
        auto alloc = Kokkos::view_alloc(label);
#else
        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);
#endif
        parts = Particles(alloc, n_reserved);

        // with possible paddings
        n_reserved = parts.stride(1);

        masks = ParticleMasks(label+"_masks", n_reserved);
        discards = ParticleMasks(label+"_discards", n_reserved);

        hparts = Kokkos::create_mirror_view(parts);
        hmasks = Kokkos::create_mirror_view(masks);
        hdiscards = Kokkos::create_mirror_view(discards);

        // id
        assign_ids(offsets_t[mpi_rank], comm);

        // valid particles
        particle_masks_initializer pmi{masks, n_valid};
        Kokkos::parallel_for(n_reserved, pmi);
    } 
}

template<>
void bunch_particles_t<double>::reserve_local(int r)
{
    if (r <= n_reserved) return;

#ifdef NO_PADDING
    Kokkos::resize(parts, r);
#else
    Kokkos::resize(Kokkos::AllowPadding, parts, r);
#endif
    n_reserved = parts.stride(1);

    Kokkos::resize(masks, n_reserved);
    Kokkos::resize(discards, n_reserved);

    hparts = Kokkos::create_mirror_view(parts);
    hmasks = Kokkos::create_mirror_view(masks);
    hdiscards = Kokkos::create_mirror_view(discards);
}

template<>
void bunch_particles_t<double>::reserve(int n, Commxx const& comm)
{
    int r = decompose_1d_local(comm, n);
    reserve_local(r);
}

template<>
void bunch_particles_t<double>::assign_ids(int train_idx, int bunch_idx)
{
    // each bunch is assined a range in the global id space
    //
    // here we assume the max number of particles that a single
    // bunch can hold is 2^32. 
    // a bunch can have 2^2, or 4 groups of particles
    // a train can have 2^16, or 65536 bunches.
    // a bunch simulator can have 2^2, or 4 trains
    //
    // so the overall bits the ids take is 32+2+16+2 = 52, less
    // than 53 which is the max integer number that can be 
    // represented in double

    int64_t base = 0;
    base |= (int64_t)train_idx << 50;
    base |= (int64_t)bunch_idx << 34;
    base |= (int64_t)group << 32;

    particle_id_assigner pia{parts, base+poffset};
    Kokkos::parallel_for(n_active, pia);
}

template<>
void bunch_particles_t<double>::inject( 
        bunch_particles_t const& o,
        karray1d_dev const& ref_st_diff,
        karray1d_dev const& tgt_st,
        karray1d_dev const& inj_st,
        double pdiff )
{
    // no nothing for empty inj bunch
    if (!o.n_active) return;

    // expand the particle array
    reserve_local(n_active + o.n_active);

    // inject
    particle_injector pi{ 
        parts, masks,
        o.parts, o.masks,
        ref_st_diff, tgt_st, inj_st, 
        n_active, pdiff
    };

    Kokkos::parallel_for(o.n_active, pi);

    // update number of particles
    n_active += o.n_active;
    n_valid += o.n_valid;
}

template<>
void
bunch_particles_t<double>::convert_to_fixed_t_lab(double p_ref, double beta)
{
    fixed_z_to_t_converter alg{parts, masks, p_ref, beta};
    Kokkos::parallel_for(n_active, alg);
}

template<>
void
bunch_particles_t<double>::convert_to_fixed_z_lab(double p_ref, double beta)
{
    fixed_t_to_z_converter alg{parts, masks, p_ref, beta};
    Kokkos::parallel_for(n_active, alg);
}

template<>
std::pair<karray2d_row, HostParticleMasks> 
bunch_particles_t<double>::get_particles_in_range(int idx, int n) const
{
    // index out of range
    if (idx == particle_index_null || idx < 0 || idx+n > n_active)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    karray2d_row_dev p("sub_p", n, 7);
    ParticleMasks pm("masks", n);

    particle_copier_many pc{parts, p, masks, pm, idx};
    Kokkos::parallel_for(n, pc);

    karray2d_row hp = create_mirror_view(p);
    Kokkos::deep_copy(hp, p);

    HostParticleMasks hpm = create_mirror_view(pm);
    Kokkos::deep_copy(hpm, pm);

    return std::make_pair(hp, hpm);
}

template<>
std::pair<karray1d_row, bool>              
bunch_particles_t<double>::get_particle(int idx) const
{
    // index out of range
    if (idx == particle_index_null || idx < 0 || idx > n_active)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    karray1d_row_dev p("particle", 7);
    ParticleMasks pm("mask", 1);

    particle_copier_one pc{parts, p, masks, pm, idx};
    Kokkos::parallel_for(1, pc);

    karray1d_row hp = create_mirror_view(p);
    Kokkos::deep_copy(hp, p);

    HostParticleMasks hpm = create_mirror_view(pm);
    Kokkos::deep_copy(hpm, pm);

    return std::make_pair(hp, hpm(0));
}

template<>
int
bunch_particles_t<double>::search_particle(int pid, int last_idx) const
{
    if (last_idx != particle_index_null)
    {
        int match = 0;
        particle_id_checker pic{parts, last_idx, pid};
        Kokkos::parallel_reduce(1, pic, match);

        if (match) return last_idx;
    }

    int idx = particle_index_null;
    particle_finder pf{parts, pid};
    Kokkos::parallel_reduce(n_active, pf, idx);

    return idx;
}

template<>
karray2d_row
bunch_particles_t<double>::get_particles_last_discarded() const
{ 
    karray2d_row_dev discarded("discarded", n_last_discarded, 7);
    karray2d_row hdiscarded = Kokkos::create_mirror_view(discarded);

    discarded_particle_mover dpm(parts, discards, discarded);
    Kokkos::parallel_for(n_active, dpm);

    Kokkos::deep_copy(hdiscarded, discarded);
    return hdiscarded;
}


#if 0
void
bunch_particles_t::set_local_num(int n)
{
    num = n;

#if 0
    // make sure the new local_num is never less than 0
    if (n < 0) n = 0;

    // re-allocate depending on the new size
    if (n <= num_aligned)
    {
        // previous values
        // int prev_local_num = num;

        // no need to resize the array, only move the pointers
        num = n;

        // update num_aligned
        num_aligned = calculate_aligned_pos(num, alignment);

        // clear the particle data (from local_num to prev_local_num)
        // note this only happens when the new local_num is smaller than the old one
        // particle_zeroer pz { parts, num };
        // Kokkos::parallel_for(prev_num - num, pz);
    }
    else
    {
        // expand the local particle array, no additional lost particle slots
        expand_local_num(num, 0);
    }
#endif
}

void
bunch_particles_t::expand_local_num(int num, int added_lost)
{
#if 0
    // keep the previous values
    int prev_local_num = local_num;
    int prev_local_num_padded = local_num_padded;

    int local_num_lost = local_num_slots - local_num_padded;
    int total_num_lost = local_num_lost + added_lost;

    double * prev_storage = storage;
    MArray2d_ref * prev_local_particles = local_particles;

    // update the pointers
    local_num = num;
    local_num_aligned = calculate_aligned_pos(local_num);
    local_num_padded = local_num + calculate_padding_size(local_num + total_num_lost);
    local_num_slots = local_num_padded + total_num_lost;

    // allocate for new storage
    storage = (double*)boost::alignment::aligned_alloc(8 * sizeof(double), local_num_slots * 7 * sizeof(double));
    local_particles = new MArray2d_ref(storage, boost::extents[local_num_slots][7], boost::fortran_storage_order());

    // copy the particle data over
    for (int i = 0; i < prev_local_num; ++i) 
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[i][j] = (*prev_local_particles)[i][j];
        }
    }

    // set the coordinates of extended and padding particles to 0
    // TODO: what should be the id for the extended particles
    for (int i = prev_local_num; i < local_num_padded; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[i][j] = 0.0;
        }
    }

    // copy over lost particles
    for (int i=0; i<local_num_lost; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[local_num_padded + i][j] = 
                (*prev_local_particles)[prev_local_num_padded + i][j];
        }
    }

    // set additional lost particle data to 0
    for (int i = local_num_padded + local_num_lost; i < local_num_slots; ++i)
    {
        for (int j=0; j<7; ++j) 
        {
            (*local_particles)[i][j] = 0.0;
        }
    }

    if (prev_storage) boost::alignment::aligned_free(prev_storage);
    if (prev_local_particles) delete prev_local_particles;
#endif
}
#endif

template<>
int
bunch_particles_t<double>::update_valid_num()
{
    int old_valid_num = n_valid;
    mask_reducer mr{masks};
    Kokkos::parallel_reduce(n_active, mr, n_valid);
    return old_valid_num;
}

template<>
int
bunch_particles_t<double>::update_total_num(Commxx const& comm)
{
    int old_total_num = n_total;
    MPI_Allreduce(&n_valid, &n_total, 1, MPI_INT, MPI_SUM, comm);
    return old_total_num;
}

template<>
void 
bunch_particles_t<double>::check_pz2_positive()
{
    checkout_particles();

    for (int p = 0; p < n_active; ++p) 
    {
        if (!hmasks(p)) continue;

        double pzop2 = (1. + hparts(p, 5)) * (1. + hparts(p, 5))
            - hparts(p, 1) * hparts(p, 1)
            - hparts(p, 3) * hparts(p, 3);

        if ( pzop2 < 0.0 )  
        {
            std::cout << "pzop^2 = " << pzop2 << std::endl;
            throw std::runtime_error( " check pz2:  pz square cannot be negative!");
        }
    }
}

template<>
void
bunch_particles_t<double>::read_file_legacy(Hdf5_file const& file, Commxx const& comm)
{
    auto dims = file.get_dims("particles");
    if (dims.size() != 2 || dims[1] != 7)
    {
        throw std::runtime_error(
                "BunchParticle::read_file_legacy(): wrong data dimensions in file");
    }

    int file_total = dims[0];
    int file_num = decompose_1d_local(comm, file_total);

    // size check
    if (n_active != file_num)
    {
        throw std::runtime_error( 
                " local_num incompatibility when initializing the bunch");
    }

    // read
    auto read_particles = file.read<karray2d_row>("particles", n_active);

    // transpose: read_particles is row major, hparts is col major
    for (int part = 0; part < n_active; ++part) 
        for (int i = 0; i < 7; ++i) 
            hparts(part, i) = read_particles(part, i);

    // check in to device mem
    checkin_particles();

    // all particles in the legacy file are valid
    // do the init after checkin because masks_initializer is performed
    // on the device memory
    n_valid = n_active;
    particle_masks_initializer pmi{masks, n_valid};
    Kokkos::parallel_for(n_reserved, pmi);
}

template<>
void
bunch_particles_t<double>::read_file(Hdf5_file const& file, Commxx const& comm)
{
    auto dims = file.get_dims(label);
    if (dims.size() != 2 || dims[1] != 7)
    {
        throw std::runtime_error(
                "BunchParticle::read_file(): wrong data dimensions in file");
    }

    int file_total = dims[0];
    int file_num = decompose_1d_local(comm, file_total);

    // allocation if capacity is smaller
    if (n_reserved < file_num)
    {
#ifdef NO_PADDING
        auto alloc = Kokkos::view_alloc(label);
#else
        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);
#endif
        parts = Particles(alloc, file_num);
        n_reserved = parts.stride(1);

        masks = ParticleMasks(label+"_masks", n_reserved);
        discards = ParticleMasks(label+"_discards", n_reserved);

        hparts = Kokkos::create_mirror_view(parts);
        hmasks = Kokkos::create_mirror_view(masks);
        hdiscards = Kokkos::create_mirror_view(discards);
    }

    // reset the pointers
    n_valid  = 0;
    n_total  = 0;
    n_active = file_num;

    // read from file
    auto read_particles = file.read<karray2d_row>(label, file_num);
    auto read_masks = file.read<HostParticleMasks>(label+"_masks", file_num);

    // transpose: read_particles is row major, hparts is col major
    for (int part = 0; part < file_num; ++part) 
    {
        for (int i = 0; i < 7; ++i) 
            hparts(part, i) = read_particles(part, i);

        hmasks(part) = read_masks(part);
        if (hmasks(part)) ++n_valid;
    }

    // now we have the actual local num, update the total number
    update_total_num(comm);

    // check in to device mem
    checkin_particles();
}

template<>
void
bunch_particles_t<double>::write_file(Hdf5_file const& file, 
        int num_part, int offset, Commxx const& comm) const
{
    int local_num_part = 0;
    int local_offset = 0;

    if (num_part == -1)
    {
        local_num_part = n_active;
        local_offset = 0;
    }
    else
    {
        local_num_part = decompose_1d_local(comm, num_part);
        local_offset = decompose_1d_local(comm, offset);
    }

    if (local_num_part < 0 || local_offset < 0 
            || local_num_part + local_offset > n_active)
    {
        throw std::runtime_error(
                "invalid num_part or offset for bunch_particles_t::write_file()");
    }

    auto parts = get_particles_in_range(local_offset, local_num_part);
    file.write_collective(label, parts.first);
    file.write_collective(label + "_masks", parts.second);
}

template<>
void 
bunch_particles_t<double>::print_particle(size_t idx, Logger& logger) const
{
    logger(LoggerV::DEBUG)
        << std::showpos << std::scientific
        << std::setprecision(8)
        << std::setw(12) << hparts(idx, 0) << ", "
        << std::setw(12) << hparts(idx, 1) << ", "
        << std::setw(12) << hparts(idx, 2) << ", "
        << std::setw(12) << hparts(idx, 3) << ", "
        << std::setw(12) << hparts(idx, 4) << ", "
        << std::setw(12) << hparts(idx, 5) << "\n"
        << std::defaultfloat << std::noshowpos
        ;
}

template<>
void
bunch_particles_t<double>::save_checkpoint_particles(Hdf5_file & file, int idx) const
{
    checkout_particles();

    std::stringstream ss;
    ss << "bunch_particles_" << label << "_parts_" << idx;
    file.write(ss.str(), hparts.data(), hparts.span(), true);

    ss.str("");
    ss << "bunch_particles_" << label << "_masks_" << idx;
    file.write(ss.str(), hmasks.data(), hmasks.span(), true);
}

template<>
void
bunch_particles_t<double>::load_checkpoint_particles(Hdf5_file & file, int idx)
{
    std::stringstream ss;
    ss << "bunch_particles_" << label << "_parts_" << idx;
    file.read(ss.str(), hparts.data(), hparts.span());

    ss.str("");
    ss << "bunch_particles_" << label << "_masks_" << idx;
    file.read(ss.str(), hmasks.data(), hmasks.span());

    checkin_particles();
}




