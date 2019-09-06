#ifndef BUNCH_PARTICLES_H
#define BUNCH_PARTICLES_H


#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"


typedef Kokkos::View<double*[7], Kokkos::LayoutLeft> Particles;
typedef Kokkos::View<const double*[7], Kokkos::LayoutLeft> ConstParticles;

typedef Particles::HostMirror HostParticles;
typedef ConstParticles::HostMirror ConstHostParticles;

class BunchParticles
{
public:

    constexpr const static int alignment = 4;
    constexpr static const int particle_index_null = -1;

private:

    /*
     * Local Particle Array Memory Layout:
     *
     *   P: particle, O: padding, L: lost particle
     *
     *   +=====+
     *   |  P  |
     *   +-----+
     *   |  P  |
     *   +-----+
     *   | ... |
     *   +=====+  <- local_num
     *   |  O  |
     *   +-----+  <- local_num_aligned
     *   |  O  |
     *   +-----+
     *   | ... |
     *   +=====+  <- local_num_padded
     *   |  L  |
     *   +-----+
     *   |  L  |
     *   +-----+
     *   | ... |
     *   +=====+  <- local_num_slots
     *
     *   * number of padding slots  = local_num_padded - local_num
     *   * number of lost particles = local_num_slots - local_num_padded
     *
     *   At bunch construction the size of padding (num_padding) is decided 
     *   such that the local_num_slots is always aligned (depending on the 
     *   vector specification, e.g., SSE or AVX or AVX512). 
     *
     *   local_num_aligned is initialized in the range [local_num, 
     *   local_num_paded], and gets adjusted everytime the local_num 
     *   changes such tht local_num_aligned is always aligned.
     *
     */

    int num;
    int num_aligned;
    int padding;
    int slots;
    int total;

    // number of particles discarded from the most recent aperture apply.
    // so the discarded particle array can be retrieved by calling:
    // get_particles_in_range(local_num_padded(), last_discarded);
    int last_discarded_;

    Commxx comm;

public:

    Particles parts;
    HostParticles hparts;

    BunchParticles(int total_num, Commxx const& comm);

    int local_num()         const { return num; }
    int local_num_aligned() const { return num_aligned; }
    int local_num_padding() const { return padding; }
    int local_num_padded()  const { return num + padding; }
    int local_num_slots()   const { return slots; }
    int total_num()         const { return total; }
    int last_discarded()    const { return last_discarded_; }

    void checkin_particles()  { Kokkos::deep_copy(parts, hparts); }
    void checkout_particles() { Kokkos::deep_copy(hparts, parts); }

    void assign_ids(int local_offset);

    void set_local_num(int num);
    void set_total_num(int num);
    void expand_local_num(int num, int added_lost);

    // update total num across the ranks and returns the old total number
    int update_total_num();

    // apply aperture operation
    template<typename AP>
    int apply_aperture(AP const& ap);

    int search_particle(int pid, int last_idx) const;
    karray1d_row get_particle(int idx) const;
    karray2d_row get_particles_in_range(int idx, int num) const;

    karray2d_row get_particles_last_discarded() const
    { return get_particles_in_range(local_num_padded(), last_discarded()); }

    void check_pz2_positive();
    void print_particle(size_t idx, Logger & logger) const;
};



namespace bunch_particles_impl
{
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
    k1d_byte discard(Kokkos::ViewAllocateWithoutInitializing("discard"), num);

    // go through each particle to see which one is been filtered out
    discard_checker<AP> dc{ap, parts, discard};
    Kokkos::parallel_reduce(num, dc, ndiscarded);

    if (ndiscarded == 0) return ndiscarded;
    std::cout << "      discarded = " << ndiscarded << "\n";

    // move discarded particles to the tail of the array
    particle_mover pm {num, padding, ndiscarded, parts, discard};
    Kokkos::parallel_for(1, pm);

    // adjust the array pointers
    last_discarded_ = ndiscarded;
    set_local_num(num - ndiscarded);

    return ndiscarded;

#if 0
    // diagnostics_loss update and write
    auto diag_loss = bunch.get_diag_loss_aperture();
    if (diag_loss)
    {
        int start_idx = bunch.get_local_num_padded() - ndiscarded;
        auto discarded = bunch.get_particles_in_range(start_idx, ndiscarded);
        //diag_loss->update(discarded);
        //diag_loss->write();
    }
#endif
}



#endif
