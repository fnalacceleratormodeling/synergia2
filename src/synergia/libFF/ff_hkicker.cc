#include "ff_hkicker.h"
#include "ff_algorithm.h"
#include "synergia/utils/simple_timer.h"

#include "ff_patterned_propagator.h"

#if 0
namespace
{
    struct thin_kicker
    {
        Particles p;
        ParticleMasks m;

        double b_hk, b_vk;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if(m(i)) 
            {
                FF_algorithm::thin_kicker_unit(p(i, 1), b_hk); 
                FF_algorithm::thin_kicker_unit(p(i, 3), b_vk); 
            }
        }
    };

    struct simple_kicker
    {
        Particles p;
        ParticleMasks m;

        double len, pref, mass, ref_cdt, b_hk, b_vk;

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

                FF_algorithm::thin_kicker_unit(p(i, 1), b_hk); 
                FF_algorithm::thin_kicker_unit(p(i, 3), b_vk); 

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
        double const* step_str;
        int steps;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if(m(i)) 
            {
                FF_algorithm::yoshida6<double, 
                    FF_algorithm::thin_kicker_unit<double>, 1>(
                            p(i, 0), p(i, 1), p(i, 2),
                            p(i, 3), p(i, 4), p(i, 5),
                            pref, mass, step_ref_cdt,
                            step_len, step_str, steps );
            }
        }
    };


    void apply_thin_kick(Bunch& bunch, ParticleGroup pg, 
            double b_hk, double b_vk)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        thin_kicker tk{parts, masks, b_hk, b_vk};
        Kokkos::parallel_for(bunch.size(pg), tk);
    }

    void apply_simple_kick(Bunch& bunch, ParticleGroup pg,
            double len, double pref, double mass, 
            double ref_cdt, double b_hk, double b_vk)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        simple_kicker sk{ parts, masks, 
            len, pref, mass, ref_cdt, b_hk, b_vk
        };

        Kokkos::parallel_for(bunch.size(pg), sk);
    }

    void apply_yoshida_kick(Bunch& bunch, ParticleGroup pg,
            double pref, double mass, double step_ref_cdt,
            double step_len, double const* step_str, int steps)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        yoshida_kicker yk{ parts, masks,
            pref, mass, step_ref_cdt, 
            step_len, step_str, steps 
        };

        Kokkos::parallel_for(bunch.size(pg), yk);
    }

    // for zero-length
    void get_reference_cdt_zero(Reference_particle& ref, 
            double hk, double vk)
    {
        double xp = ref.get_state()[Bunch::xp];
        double yp = ref.get_state()[Bunch::yp];

        FF_algorithm::thin_kicker_unit(xp, hk);
        FF_algorithm::thin_kicker_unit(yp, vk);

        ref.set_state_xp(xp);
        ref.set_state_yp(yp);
        ref.set_state_cdt(0.0);
    }

    // non zero length
    double get_reference_cdt_simple(Reference_particle& ref, 
            double len, double hk, double vk)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();

        FF_algorithm::drift_unit(
                st[0], st[1], st[2], 
                st[3], cdt, st[5], 
                len* 0.5, pref, mass, 0.0);

        // for >0 length, hk,vk is the strength/length of the kick
        FF_algorithm::thin_kicker_unit(st[1], hk * len);
        FF_algorithm::thin_kicker_unit(st[3], vk * len);

        FF_algorithm::drift_unit(
                st[0], st[1], st[2], 
                st[3], cdt, st[5], 
                len* 0.5, pref, mass, 0.0);

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }

    // non zero length
    double get_reference_cdt_yoshida(Reference_particle& ref, 
            double len, double hk, double vk, int steps)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();

        // steps comes from base class, set in apply method
        double step_len = len / steps;

        // for >0 length, hk,vk is the strength/length of the kick
        double step_str[2] = { hk*step_len, vk*step_len };

        // propagate
        FF_algorithm::yoshida6<double, 
            FF_algorithm::thin_kicker_unit<double>, 1>( 
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5],
                    pref, mass, 0.0,
                    step_len, step_str, steps );

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }
}
#endif

void FF_hkicker::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error(
            "libFF::hkicker for JetParticles has yet to be implemented");
}

void FF_hkicker::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    scoped_simple_timer timer("libFF_hkicker");

    auto const& elem = slice.get_lattice_element();
    const double length = slice.get_right() - slice.get_left();

    // hk and vk are the hk/vk under lattice reference momentum
    double    l = elem.get_double_attribute("l");
    double  hk0 = elem.get_double_attribute("hkick");
    double  vk0 = elem.get_double_attribute("vkick");
    double tilt = elem.get_double_attribute("tilt");

    // tilt
    double hk = cos(tilt)*hk0 - sin(tilt)*vk0;
    double vk = sin(tilt)*hk0 + cos(tilt)*vk0;

    //double k = elem.get_double_attribute(k_attr, 0.0);

    auto& ref_lattice = bunch.get_design_reference_particle();
    auto const& ref_bunch = bunch.get_reference_particle();

    double plattice = ref_lattice.get_momentum();
    double pbunch = ref_bunch.get_momentum();

    // scale is to scale the kick strength defined relative to the lattice momentum to
    // the scale of the bunch particles defined relative to the bunch momentum
    double scale = plattice/pbunch;

    // kick strength is defined as momentum change/reference momentum
    double b_hk = hk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;
    double b_vk = vk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;

    double k[2] = {hk, vk};
    double sk[2] = {b_hk, b_vk};

    using pp = FF_patterned_propagator<double, 1,
          FF_algorithm::thin_kicker_unit<double>>;

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference 
        // particle which is relative to the p-lattice.
        // also update the reference particle
        pp::get_reference_cdt_zero(ref_lattice, k);

        pp::apply_thin_kick(bunch, ParticleGroup::regular, sk);
        pp::apply_thin_kick(bunch, ParticleGroup::spectator, sk);
    }
    else
    {
        // simple drift-kick-drift scheme
        double simple_d = slice
            .get_lattice_element()
            .get_double_attribute("simple", 0.0);

        bool simple = fabs(simple_d) > 1e-16;

        // strength per unit length
#if 0
        double hk_pul = hk / l;
        double vk_pul = vk / l;
#endif

        k[0] = k[0]/l;
        k[1] = k[1]/l;

        sk[0] = sk[0]/l;
        sk[1] = sk[1]/l;

        double pref = bunch.get_reference_particle().get_momentum();
        double mass = bunch.get_mass();

        if (simple)
        {
            // use un-scaled k_pul
            double ref_cdt = pp::get_reference_cdt_simple(
                    ref_lattice, length, k);

            pp::apply_simple_kick(bunch, ParticleGroup::regular, 
                    pref, mass, ref_cdt, length, sk);

            pp::apply_simple_kick(bunch, ParticleGroup::spectator, 
                    pref, mass, ref_cdt, length, sk);
        }
        else
        {
            // yoshida steps
            steps = (int)slice
                .get_lattice_element()
                .get_double_attribute("yoshida_steps", 4.0);

            // use un-scaled k_pul
            double ref_cdt = pp::get_reference_cdt_yoshida(
                    ref_lattice, length, k, steps, true);

#if 0
            double b_hk_pul = b_hk / l;
            double b_vk_pul = b_vk / l;

            double step_ref_cdt = ref_cdt / steps;
            double step_len = length / steps;
            double step_str[2] = { b_hk_pul * step_len, b_vk_pul * step_len };
#endif

            pp::apply_yoshida_kick(bunch, ParticleGroup::regular,
                    pref, mass, ref_cdt, length, sk, steps);

            pp::apply_yoshida_kick(bunch, ParticleGroup::spectator,
                    pref, mass, ref_cdt, length, sk, steps);
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
}

