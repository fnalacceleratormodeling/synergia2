#ifndef BUNCH_PARTICLES_H
#define BUNCH_PARTICLES_H


#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/gsvector.h"

#include <cereal/cereal.hpp>

enum class ParticleGroup
{
    regular = 0,
    spectator = 1
};

typedef Kokkos::View<double*[7], Kokkos::LayoutLeft> Particles;
typedef Kokkos::View<const double*[7], Kokkos::LayoutLeft> ConstParticles;

typedef Particles::HostMirror HostParticles;
typedef ConstParticles::HostMirror ConstHostParticles;

typedef Kokkos::View<uint8_t*> ParticleMasks;
typedef Kokkos::View<const uint8_t*> ConstParticleMasks;

typedef ParticleMasks::HostMirror HostParticleMasks;
typedef ConstParticleMasks::HostMirror ConstHostParticleMasks;

// serialization
namespace cereal
{
    // particles
    template<class AR>
    void save(AR & ar, Particles const& p)
    {
        ar(p.label(), p.stride(1));
    }

    template<class AR>
    void load(AR & ar, Particles & p)
    {
        std::string label;
        int slots;
        ar(label, slots);

        auto alloc = Kokkos::view_alloc(label, Kokkos::AllowPadding);
        p = Particles(alloc, slots);

        if (p.stride(1) != slots)
            throw std::runtime_error("inconsistent padding while loading particles");
    }

    // masks
    template<class AR>
    void save(AR & ar, ParticleMasks const& p)
    {
        ar(p.label(), p.extent(0));
    }

    template<class AR>
    void load(AR & ar, ParticleMasks & p)
    {
        std::string label;
        int slots;

        ar(label, slots);
        p = ParticleMasks(label, slots);
    }
}


template<class PART>
class bunch_particles_t
{
    using PG = ParticleGroup;

public:

    constexpr static const int particle_index_null = -1;

    using part_t = PART;

    using parts_t = typename Kokkos::View<PART*[7], Kokkos::LayoutLeft>;
    using const_parts_t = typename Kokkos::View<const PART*[7], Kokkos::LayoutLeft>;

    using host_parts_t = typename parts_t::HostMirror;
    using const_host_parts_t = typename const_parts_t::HostMirror;

    using gsv_t = typename std::conditional<
        std::is_floating_point<PART>::value, 
        GSVector, Vec<PART>>::type;

private:

    /*
     * Local Particle Array Memory Layout:
     *
     *   P: regular particle, 
     *   I: invalid (lost) particle
     *   R: reserved particle
     *
     *    part       mask
     *   +=====+    +=====+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  I  |    |  0  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  I  |    |  0  |
     *   +-----+    +-----+
     *   |  I  |    |  0  |
     *   +-----+    +-----+
     *   |  P  |    |  1  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +-----+    +-----+
     *   |  R  |    |  0  |
     *   +=====+    +=====+ 
     *
     *   num_valid    = 5   -- num_valid()
     *   num_active   = 8   -- num_active() / size()
     *   num_reserved = 12  -- num_reserved() / capacity()
     *
     *   when allocating an particle array, the padding flag of the 
     *   Kokkos::View object is turned on. So the 'num_reserved' is 
     *   the actual array size in the first dimension.
     *
     *   'num_valid' is the number of remainig local particles in
     *   the bunch
     *
     *   'num_active' is the number of local particles including lost
     *   ones. This one should be used when looping through the 
     *   particles array
     *
     */

    // particle group (regular or spectator)
    ParticleGroup group;
    std::string label;

    // see the memory layout
    int n_valid;
    int n_active;
    int n_reserved;

    // sum of n_valid on all ranks
    int n_total; 

    // number of particles discarded from the most recent aperture apply.
    int n_last_discarded;

    // particle offset in cases where a single bunch span across
    // multiple ranks
    int poffset;


public:

    parts_t parts;
    ParticleMasks masks;
    ParticleMasks discards;

    host_parts_t hparts;
    HostParticleMasks hmasks;
    HostParticleMasks hdiscards;

public:

    // particles array will be allocated to the size of reserved_num
    // if reserved is less than total, it will be set to the total_num
    bunch_particles_t( 
            ParticleGroup pg,
            int total_num, 
            int reserved_num, 
            Commxx const& comm );

    // constructor for trigon particles
    template<typename U = PART>
    bunch_particles_t(
            int total_num,
            Commxx const& comm, 
            typename std::enable_if<U::is_trigon>::type* = 0);

    int num_valid()          const { return n_valid; }
    int num_active()         const { return n_active; }
    int num_reserved()       const { return n_reserved; }
    int num_total()          const { return n_total; }
    int num_last_discarded() const { return n_last_discarded; }

    // getters with names more consistent with std containers
    int size()     const { return n_active; }
    int capacity() const { return n_reserved; }

    // copy particles/masks between host and device memories
    void checkin_particles() const
    { Kokkos::deep_copy(parts, hparts); Kokkos::deep_copy(masks, hmasks); }

    void checkout_particles() const
    { Kokkos::deep_copy(hparts, parts); Kokkos::deep_copy(hmasks, masks); }

    // change capacity (can only increase)
    // n is the new cap across the entire bunch
    void reserve(int n, Commxx const& comm);

    // change local capacity (can only increase)
    // n is the new local capacity
    void reserve_local(int n);

    // inject with 
    void inject(bunch_particles_t const& o,
            karray1d_dev const& ref_st_diff,
            karray1d_dev const& tgt_st,
            karray1d_dev const& inj_st,
            double pdiff );

    // convert between fixed z lab and fixed t lab
    void convert_to_fixed_t_lab(double p_ref, double beta);
    void convert_to_fixed_z_lab(double p_ref, double beta);

#if 0
    void set_total_num(int num);
    void expand_local_num(int num, int added_lost);
#endif

    // update the valid num from the masks and return the old valid num
    int update_valid_num();

    // update total num across the ranks and returns the old total number
    int update_total_num(Commxx const& comm);

    // apply aperture operation
    template<typename AP>
    int apply_aperture(AP const& ap);

    // search/get particle(s)
    int search_particle(int pid, int last_idx) const;

    std::pair<karray1d_row, bool>              
    get_particle(int idx) const;

    std::pair<karray2d_row, HostParticleMasks> 
    get_particles_in_range(int idx, int num) const;

    karray2d_row 
    get_particles_last_discarded() const;

    void check_pz2_positive();
    void print_particle(size_t idx, Logger& logger) const;

    // read from a hdf5 file. total_num of current bunch must be the same 
    // as the one stored in the particle file
    void read_file_legacy(Hdf5_file const& file, Commxx const& comm);

    void read_file (Hdf5_file const& file, Commxx const& comm);
    void write_file(Hdf5_file const& file, int num_part, int offset, Commxx const& comm) const;

    // checkpoint save/load
    void save_checkpoint_particles(Hdf5_file & file, int idx) const;
    void load_checkpoint_particles(Hdf5_file & file, int idx);

    // assign ids cooperatively
    void assign_ids(int train_idx, int bunch_idx);

private:

    void assign_ids(int local_offset, Commxx const& comm);

    // serialization
    friend class cereal::access;

    template<class AR>
    void save(AR & ar) const
    { 
        ar(CEREAL_NVP(label));
        ar(CEREAL_NVP(n_valid));
        ar(CEREAL_NVP(n_active));
        ar(CEREAL_NVP(n_reserved));
        ar(CEREAL_NVP(n_total));

        ar(CEREAL_NVP(parts));
        ar(CEREAL_NVP(masks));
        ar(CEREAL_NVP(discards));
    }

    template<class AR>
    void load(AR & ar)
    { 
        ar(CEREAL_NVP(label));
        ar(CEREAL_NVP(n_valid));
        ar(CEREAL_NVP(n_active));
        ar(CEREAL_NVP(n_reserved));
        ar(CEREAL_NVP(n_total));

        ar(CEREAL_NVP(parts));
        ar(CEREAL_NVP(masks));
        ar(CEREAL_NVP(discards));

        // construct the host array first
        hparts = Kokkos::create_mirror_view(parts);
        hmasks = Kokkos::create_mirror_view(masks);
        hdiscards = Kokkos::create_mirror_view(discards);
    }
};


// instantiated with a double is still the default BunchParticles
typedef bunch_particles_t<double> BunchParticles;


// implementations
namespace bunch_particles_impl
{
    template<class AP>
    struct discard_applier
    {
        typedef int value_type;

        AP ap;
        ConstParticles parts;
        ParticleMasks masks;
        ParticleMasks discards;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& discarded) const
        {
            discards(i) = 0;

            if (masks(i) && ap.discard(parts, i))
            {
                discards(i) = 1;
                masks(i) = 0;
                ++discarded;
            }
        }
    };
}


template<>
template<typename AP>
inline int bunch_particles_t<double>::apply_aperture(AP const& ap)
{
    using namespace bunch_particles_impl;

    int ndiscarded = 0;

    // go through each particle to see which one is been filtered out
    discard_applier<AP> da{ap, parts, masks, discards};
    Kokkos::parallel_reduce(n_active, da, ndiscarded);

    //std::cout << "      discarded = " << ndiscarded << "\n";
    n_last_discarded = ndiscarded;
    n_valid -= ndiscarded;

    return ndiscarded;
}


template<typename PART>
template<typename U>
inline bunch_particles_t<PART>::bunch_particles_t(
        int total_num,
        Commxx const& comm,
        typename std::enable_if<U::is_trigon>::type*)
    : group(PG::regular)
    , label("trigons")
    , n_valid((!comm.is_null() && comm.rank() == 0) ? total_num : 0)
    , n_active(n_valid)
    , n_reserved(n_valid)
    , n_total(n_valid)
    , n_last_discarded(0)
    , poffset(0)
    , parts(parts_t("trigon", n_reserved))
    , masks(ParticleMasks("trigon_masks", n_reserved))
    , discards()
    , hparts(Kokkos::create_mirror_view(parts))
    , hmasks(Kokkos::create_mirror_view(masks))
    , hdiscards()
{
    for(int i=0; i<n_valid; ++i) hmasks(i) = 1;
    Kokkos::deep_copy(masks, hmasks);
}

#endif
