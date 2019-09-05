
#include <iomanip>

#include "synergia/bunch/bunch_particles.h"
#include "synergia/utils/parallel_utils.h"



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
        int offset;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { parts(i, 6) = i + offset; }
    };


    /// For a given number of particles, returns the next alignement position
    int calculate_aligned_pos(int num, int alignment)
    {
        if (alignment <= 0)
            throw std::runtime_error("calculate_aligned_pos() invalid particle_alignment value");

        if (num < 0)
            throw std::runtime_error("calculate_aligned_pos() invalid num value");

        if (num == 0) return 0;
        if (num % alignment == 0) return num;
        return num + alignment - (num % alignment);
    }

    /// For a given number of particles, returns the needed size of padding
    /// for the particle array to be aligned
    int calculate_padding_size(int num, int alignment)
    {
        if (alignment <= 0)
            throw std::runtime_error("calculate_padding_size() invalid particle_alignment value");

        if (num < 0)
            throw std::runtime_error("calculate_padding_size() invalid num value");

        if (num == 0) return 0;
        if (num % alignment == 0) return alignment;
        return alignment * 2 - num % alignment;
    }


    // Kokkos functors
    struct particle_copier_many
    {
        ConstParticles src;
        karray2d_row_dev dst;
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
        }
    };

    struct particle_copier_one
    {
        ConstParticles src;
        karray1d_row_dev dst;
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

}


BunchParticles::BunchParticles(int total_num, Commxx const& comm)
    : num(0)
    , num_aligned(0)
    , padding(0)
    , slots(0)
    , total(total_num)
    , comm(comm)
    , parts("particles", slots)
    , hparts(Kokkos::create_mirror_view(parts))
{
    if (!comm.is_null() && total) 
    {
        std::vector<int> offsets(comm.size());
        std::vector<int> counts(comm.size());
        decompose_1d(comm, total, offsets, counts);

        num         = counts[comm.rank()];
        num_aligned = calculate_aligned_pos(num, alignment);
        padding     = calculate_padding_size(num, alignment);
        slots       = num + padding;

        // allocate
        Kokkos::resize(parts, slots);
        hparts = Kokkos::create_mirror_view(parts);

        // id
        assign_ids(offsets[comm.rank()]);
    } 
}


void BunchParticles::assign_ids(int local_offset)
{
    int request_num = (comm.rank() == 0) ? total : 0;
    int global_offset = Particle_id_offset::get(request_num, comm);

    particle_id_assigner pia{parts, local_offset + global_offset};
    Kokkos::parallel_for(num, pia);
}

karray2d_row
BunchParticles::get_particles_in_range(int idx, int n) const
{
    karray2d_row_dev p("sub_p", n, 7);
    particle_copier_many pc{parts, p, idx};
    Kokkos::parallel_for(n, pc);

    karray2d_row hp = create_mirror_view(p);
    Kokkos::deep_copy(hp, p);
    return hp;
}

karray1d_row
BunchParticles::get_particle(int idx) const
{
    // index out of range
    if (idx == particle_index_null || idx < 0 || idx > slots)
        throw std::runtime_error("Bunch::get_particle() index out of range");

    karray1d_row_dev p("particle", 7);
    particle_copier_one pc{parts, p, idx};
    Kokkos::parallel_for(1, pc);

    karray1d_row hp = create_mirror_view(p);
    Kokkos::deep_copy(hp, p);
    return hp;
}

int
BunchParticles::search_particle(int pid, int last_idx) const
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
    Kokkos::parallel_reduce(num, pf, idx);

    return idx;
}

void
BunchParticles::set_local_num(int n)
{
    // make sure the new local_num is never less than 0
    if (n < 0) n = 0;

    // re-allocate depending on the new size
    if (n <= num_aligned)
    {
        // previous values
        int prev_local_num = num;

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
}

void
BunchParticles::expand_local_num(int num, int added_lost)
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

int
BunchParticles::update_total_num()
{
    int old_total_num = total;
    MPI_Allreduce(&num, &total, 1, MPI_INT, MPI_SUM, comm);
    return old_total_num;
}

void
BunchParticles::set_total_num(int totalnum)
{
#if 0
    int old_total_num = total_num;
    total_num = totalnum;

    if (old_total_num != 0) 
    {
        real_num = (total_num * real_num) / old_total_num;
    } 
    else 
    {
        real_num = 0.0;
    }
#endif
}

void 
BunchParticles::check_pz2_positive()
{
    checkout_particles();

    for (int p = 0; p < num; ++p) 
    {
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

void 
BunchParticles::print_particle(size_t idx, Logger & logger) const
{
    auto flags = std::cout.flags();

    logger(LoggerV::DEBUG)
        << std::resetiosflags(std::ios::fixed)
        << std::setiosflags(std::ios::showpos | std::ios::scientific)
        << std::setprecision(8)
        << std::setw(12) << hparts(idx, 0) << ", "
        << std::setw(12) << hparts(idx, 1) << ", "
        << std::setw(12) << hparts(idx, 2) << ", "
        << std::setw(12) << hparts(idx, 3) << ", "
        << std::setw(12) << hparts(idx, 4) << ", "
        << std::setw(12) << hparts(idx, 5) << "\n"
        ;

    std::cout.flags(flags);
}




