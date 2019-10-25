#ifndef BUNCH_PARTICLES_H
#define BUNCH_PARTICLES_H


#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"

typedef Kokkos::View<double*[7], Kokkos::LayoutLeft> Particles;
typedef Kokkos::View<const double*[7], Kokkos::LayoutLeft> ConstParticles;

typedef Particles::HostMirror HostParticles;
typedef ConstParticles::HostMirror ConstHostParticles;

typedef Kokkos::View<uint8_t*> ParticleMasks;
typedef Kokkos::View<const uint8_t*> ConstParticleMasks;

typedef ParticleMasks::HostMirror HostParticleMasks;
typedef ConstParticleMasks::HostMirror ConstHostParticleMasks;


class BunchParticles
{
public:

    constexpr const static int alignment = 4;
    constexpr static const int particle_index_null = -1;

private:

    /*
     * Local Particle Array Memory Layout:
     *
     *   P: particle, I: invalid particle
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
     *   +=====+    +=====+ 
     *
     *   slots = 8   (num of slots)
     *   num = 5     (num of local particles)
     *
     *   when allocating an particle array, the padding flag of the Kokkos::View
     *   object is turned on. So the 'slots' is the actual array size in the 
     *   first dimension.
     *
     *   'num' is the actual number of local particles.
     *
     */


    int num;
    int slots;
    int total;

    // number of particles discarded from the most recent aperture apply.
    int last_discarded_;

    Commxx comm;

public:

    Particles parts;
    ParticleMasks masks;

    HostParticles hparts;
    HostParticleMasks hmasks;

    BunchParticles(int total_num, Commxx const& comm);

    int local_num()         const { return num; }
    int local_num_slots()   const { return slots; }
    int total_num()         const { return total; }
    int last_discarded()    const { return last_discarded_; }

    void checkin_particles()  { Kokkos::deep_copy(parts, hparts); }
    void checkout_particles() { Kokkos::deep_copy(hparts, parts); }

    void assign_ids(int local_offset);

    void set_local_num(int num);
    void set_total_num(int num);
    void expand_local_num(int num, int added_lost);

    // read from a hdf5 file. total_num of current bunch must be the same 
    // as the one stored in the particle file
    void read_file(std::string const& filename);

    // update total num across the ranks and returns the old total number
    int update_total_num();

    // apply aperture operation
    template<typename AP>
    int apply_aperture(AP const& ap);

    int search_particle(int pid, int last_idx) const;
    karray1d_row get_particle(int idx) const;
    karray2d_row get_particles_in_range(int idx, int num) const;
    karray2d_row get_particles_last_discarded() const;

    void check_pz2_positive();
    void print_particle(size_t idx, Logger & logger) const;
};


namespace bunch_particles_impl
{
    template<class AP>
    struct discard_applier
    {
        typedef int value_type;

        AP ap;
        ConstParticles parts;
        ParticleMasks masks;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& discarded) const
        {
            if (masks(i) && ap.discard(parts, i))
            {
                masks(i) = 0;
                ++discarded;
            }
        }
    };

    using k1d_byte = Kokkos::View<uint8_t*>;

    template<class AP>
    struct discard_checker
    {
        typedef int value_type;

        AP ap;
        ConstParticles parts;
        k1d_byte discard;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i, int& discarded) const
        {
            if (ap.discard(parts, i))
            {
                discard[i] = 1;
                ++discarded;
            }
            else
            {
                discard[i] = 0;
            }
        }
    };

    struct particle_mover
    {
        int nparts;
        int padding;
        int ndiscarded;

        Particles parts;
        k1d_byte discard;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        {
            // move all the discarded particles to the tail
            int head = 0;
            int tail = nparts - 1;

            do
            {
                while (!discard[head] && head<tail) ++head;
                if (head >= tail) break;

                while ( discard[tail] && tail>head) --tail;
                if (head >= tail) break;

                for (int idx=0; idx<7; ++idx)
                {
                    double tmp = parts(head, idx);
                    parts(head, idx) = parts(tail, idx);
                    parts(tail, idx) = tmp;
                }

                ++head;
                --tail;

            } while (head < tail);

            // move some lost particles over to the padding area
            int padded = nparts + padding;
            int np = ndiscarded < padding ? ndiscarded : padding;

            for (int p=0; p<np; ++p)
            {
                // pl: position of next lost particle
                // pp: position of next padding slot
                int pl = nparts - ndiscarded + p;
                int pp = padded - 1 - p;

                for (int idx=0; idx<7; ++idx)
                {
                    // copy the lost particle over to the padding slot
                    // makes pl the new padding slot
                    parts(pp, idx) = parts(pl, idx);
                    parts(pl, idx) = 0.0;
                }
            }
        }
    };
}



template<typename AP>
inline int BunchParticles::apply_aperture(AP const& ap)
{
    using namespace bunch_particles_impl;

    int ndiscarded = 0;

    // go through each particle to see which one is been filtered out
    discard_applier<AP> da{ap, parts, masks};
    Kokkos::parallel_reduce(slots, da, ndiscarded);

    //std::cout << "      discarded = " << ndiscarded << "\n";
    set_local_num(num - ndiscarded);
    return ndiscarded;

#if 0
    int ndiscarded = 0;
    k1d_byte discard(Kokkos::ViewAllocateWithoutInitializing("discard"), num);

    // go through each particle to see which one is been filtered out
    discard_checker<AP> dc{ap, parts, discard};
    Kokkos::parallel_reduce(num, dc, ndiscarded);

    if (ndiscarded == 0) return ndiscarded;
    // std::cout << "      discarded = " << ndiscarded << "\n";

    // move discarded particles to the tail of the array
    particle_mover pm {num, padding, ndiscarded, parts, discard};
    Kokkos::parallel_for(1, pm);

    // adjust the array pointers
    last_discarded_ = ndiscarded;
    set_local_num(num - ndiscarded);

    return ndiscarded;
#endif
}



#endif
