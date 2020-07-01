#ifndef SYNERGIA_LIBFF_PATTERNED_PROPAGATOR_H
#define SYNERGIA_LIBFF_PATTERNED_PROPAGATOR_H

template<
    typename T,
    int COMP,
    void(KF)(T const&, T&, T const&, T&, T const&, double const*) >
struct FF_patterned_propagator
{
    struct thin_kicker
    {
        Particles p;
        ParticleMasks m;

        double k[COMP*2];

        thin_kicker(Particles p, ParticleMasks m, 
                double const* str)
        : p(p), m(m)
        { for(int i=0; i<COMP*2; ++i) k[i] = str[i]; }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if(m(i)) 
                KF(p(i,0), p(i,1), p(i,2), p(i,3), p(i,5), k);
        }
    };

    struct simple_kicker
    {
        Particles p;
        ParticleMasks m;

        double pref, mass, ref_cdt, len;
        double k[COMP*2];

        simple_kicker(Particles p, ParticleMasks m, 
                double pref, double mass, double ref_cdt, double len,
                double const* str)
        : p(p), m(m), pref(pref), mass(mass)
        , ref_cdt(ref_cdt), len(len)
        { for(int i=0; i<COMP*2; ++i) k[i] = str[i]; }

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if(m(i)) 
            {
                FF_algorithm::drift_unit(
                        p(i, 0), p(i, 1), p(i, 2), 
                        p(i, 3), p(i, 4), p(i, 5),
                        len*0.5, pref, mass, ref_cdt*0.5 
                );

                KF(p(i,0), p(i,1), p(i,2), p(i,3), p(i,5), k);

                FF_algorithm::drift_unit(
                        p(i, 0), p(i, 1), p(i, 2), 
                        p(i, 3), p(i, 4), p(i, 5),
                        len*0.5, pref, mass, ref_cdt*0.5 
                );
            }
        }
    };


    struct yoshida_kicker
    {
        Particles p;
        ParticleMasks m;

        double pref, mass, step_ref_cdt, step_len;
        double step_k[COMP*2];
        int steps;

        yoshida_kicker(Particles p, ParticleMasks m, 
                double pref, double mass, double step_ref_cdt, 
                double step_len, double const* step_str, int steps)
        : p(p), m(m), pref(pref), mass(mass)
        , step_ref_cdt(step_ref_cdt), step_len(step_len)
        , steps(steps)
        { for(int i=0; i<COMP*2; ++i) step_k[i] = step_str[i]; }


        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if(m(i)) 
            {
                FF_algorithm::yoshida6<double, KF, COMP>( 
                        p(i, 0), p(i, 1), p(i, 2),
                        p(i, 3), p(i, 4), p(i, 5),
                        pref, mass, step_ref_cdt,
                        step_len, step_k, steps );
            }
        }
    };

    static void apply_thin_kick(Bunch& bunch, ParticleGroup pg, 
            double const* k)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        thin_kicker tk(parts, masks, k);
        Kokkos::parallel_for(bunch.size(pg), tk);
    }

    static void apply_simple_kick(Bunch& bunch, ParticleGroup pg,
            double pref, double mass, double ref_cdt, 
            double len, double const* k)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        simple_kicker sk( parts, masks, 
            pref, mass, ref_cdt, len, k
        );

        Kokkos::parallel_for(bunch.size(pg), sk);
    }


    static void apply_yoshida_kick(Bunch& bunch, ParticleGroup pg,
            double pref, double mass, double ref_cdt,
            double len, double const* str, int steps)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        double step_len = len / steps;
        double step_ref_cdt = ref_cdt / steps;
        double step_str[2*COMP];
        for(int i=0; i<2*COMP; ++i) step_str[i] = str[i]*step_len;

        yoshida_kicker yk( parts, masks,
            pref, mass, step_ref_cdt, 
            step_len, step_str, steps 
        );

        Kokkos::parallel_for(bunch.size(pg), yk);
    }

    static void get_reference_cdt_zero(
            Reference_particle& ref, 
            double const* k)
    {
        // propagate the bunch design reference particle
        double x  = ref.get_state()[Bunch::x];
        double xp = ref.get_state()[Bunch::xp];
        double y  = ref.get_state()[Bunch::y];
        double yp = ref.get_state()[Bunch::yp];
        double dpop = ref.get_state()[Bunch::dpop];

        KF(x, xp, y, yp, dpop, k);

        ref.set_state_xp(xp);
        ref.set_state_yp(yp);
        ref.set_state_cdt(0.0);
    }

    static double get_reference_cdt_simple(
            Reference_particle& ref, 
            double len, double const* k)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();

        // for >0 length, hk,vk is the strength/length of the kick
        double str[2*COMP];
        for(int i=0; i<2*COMP; ++i) str[i] = k[i]*len;

        FF_algorithm::drift_unit(
                st[0], st[1], st[2], 
                st[3], cdt, st[5], 
                len* 0.5, pref, mass, 0.0);

        KF(st[0], st[1], st[2], st[3], st[5], str);

        FF_algorithm::drift_unit(
                st[0], st[1], st[2], 
                st[3], cdt, st[5], 
                len* 0.5, pref, mass, 0.0);

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }


    // non zero length
    static double get_reference_cdt_yoshida(
            Reference_particle& ref, 
            double len, double const* k, int steps, 
            bool use_drift = false)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();

        // steps comes from base class, set in apply method
        double step_len = len / steps;

        // for >0 length, hk,vk is the strength/length of the kick
        double step_str[2*COMP];
        for(int i=0; i<2*COMP; ++i) step_str[i] = k[i]*step_len;

        // propagate
        if (use_drift)
        {
            FF_algorithm::drift_unit(
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5], 
                    len, pref, mass, 0.0);

        }
        else
        {
            // hardcoded to use 6th order
            FF_algorithm::yoshida<double, KF, 6, COMP>( 
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5],
                    pref, mass, 0.0,
                    step_len, step_str, steps );
        }

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }

};

#endif
