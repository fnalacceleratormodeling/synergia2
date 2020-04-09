#include "ff_hkicker.h"
#include "ff_algorithm.h"
#include "synergia/utils/simple_timer.h"


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
    //double b_k = k * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;
    double b_hk = hk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;
    double b_vk = vk * ( ref_bunch.get_charge() / ref_lattice.get_charge() ) * scale;

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference 
        // particle which is relative to the p-lattice.
        // also update the reference particle
        get_reference_cdt_zero(ref_lattice, hk, vk);

        apply_thin_kick(bunch, ParticleGroup::regular, b_hk, b_vk);
        apply_thin_kick(bunch, ParticleGroup::spectator, b_hk, b_vk);
    }
    else
    {
        // simple drift-kick-drift scheme
        double simple_d = slice
            .get_lattice_element()
            .get_double_attribute("simple", 0.0);

        bool simple = fabs(simple_d) > 1e-16;

        // strength per unit length
        double hk_pul = hk / l;
        double vk_pul = vk / l;

        double pref = bunch.get_reference_particle().get_momentum();
        double mass = bunch.get_mass();

        if (simple)
        {
            // use un-scaled k_pul
            double ref_cdt = get_reference_cdt_simple(
                    ref_lattice, length, hk_pul, vk_pul);

            apply_simple_kick(bunch, ParticleGroup::regular, 
                    length, pref, mass, ref_cdt, b_hk, b_vk);

            apply_simple_kick(bunch, ParticleGroup::spectator, 
                    length, pref, mass, ref_cdt, b_hk, b_vk);
        }
        else
        {
            // yoshida steps
            steps = (int)slice
                .get_lattice_element()
                .get_double_attribute("yoshida_steps", 4.0);

            // use un-scaled k_pul
            double ref_cdt = get_reference_cdt_yoshida(
                    ref_lattice, length, hk_pul, vk_pul, steps);

            double b_hk_pul = b_hk / l;
            double b_vk_pul = b_vk / l;

            double step_ref_cdt = ref_cdt / steps;
            double step_len = length / steps;
            double step_str[2] = { b_hk_pul * step_len, b_vk_pul * step_len };

            apply_yoshida_kick(bunch, ParticleGroup::regular,
                    pref, mass, step_ref_cdt, step_len, step_str, steps);

            apply_yoshida_kick(bunch, ParticleGroup::spectator,
                    pref, mass, step_ref_cdt, step_len, step_str, steps);

        }

        bunch.get_reference_particle().increment_trajectory(length);
    }

#if 0

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    const int gsvsize = GSVector::size();

    double * RESTRICT xa;
    double * RESTRICT xpa;
    double * RESTRICT ya;
    double * RESTRICT ypa;
    double * RESTRICT cdta;
    double * RESTRICT dpopa;

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference particle which
        // relative to plattice.
        // also update the reference particle
        double reference_cdt = get_reference_cdt(0.0, k, ref_lattice, false);

        // real particles
        {
            bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

            const int num_blocks = local_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector xp(&xpa[part]);
                FF_algorithm::thin_kicker_unit(xp, b_k);
                xp.store(&xpa[part]);
            }

            for (int part = block_last; part < local_num; ++part)
            {
                double xp(xpa[part]);
                FF_algorithm::thin_kicker_unit(xp, b_k);
                xpa[part] = xp;
            }
        }

        // spectator particles
        {
            bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

            const int num_blocks = local_s_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector xp(&xpa[part]);
                FF_algorithm::thin_kicker_unit(xp, b_k);
                xp.store(&xpa[part]);
            }

            for (int part = block_last; part < local_s_num; ++part)
            {
                double xp(xpa[part]);
                FF_algorithm::thin_kicker_unit(xp, b_k);
                xpa[part] = xp;
            }
        }
    }
    else
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // simple drift-kick-drift scheme
        double simple_d = slice.get_lattice_element().get_double_attribute("simple", 0.0);
        bool simple = fabs(simple_d) > 1e-16;

        // strength per unit length
        double k_pul = k / l;

        double pref = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k_pul, ref_lattice, simple);

        if (simple)
        {
            // real particles
            {
                bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

                const int num_blocks = local_num / gsvsize;
                const int block_last = num_blocks * gsvsize;

                #pragma omp parallel for
                for (int part = 0; part < block_last; part += gsvsize)
                {
                    GSVector x(&xa[part]);
                    GSVector xp(&xpa[part]);
                    GSVector y(&ya[part]);
                    GSVector yp(&ypa[part]);
                    GSVector cdt(&cdta[part]);
                    GSVector dpop(&dpopa[part]);

                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_kicker_unit(xp, b_k);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                    x.store(&xa[part]);
                    xp.store(&xpa[part]);
                    y.store(&ya[part]);
                    yp.store(&ypa[part]);
                    cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_kicker_unit(xp, b_k);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                    xa[part] = x;
                    xpa[part] = xp;
                    ya[part] = y;
                    ypa[part] = yp;
                    cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }
 
            // spectator particles
            {
                bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

                const int num_blocks = local_s_num / gsvsize;
                const int block_last = num_blocks * gsvsize;

                #pragma omp parallel for
                for (int part = 0; part < block_last; part += gsvsize)
                {
                    GSVector x(&xa[part]);
                    GSVector xp(&xpa[part]);
                    GSVector y(&ya[part]);
                    GSVector yp(&ypa[part]);
                    GSVector cdt(&cdta[part]);
                    GSVector dpop(&dpopa[part]);

                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_kicker_unit(xp, b_k);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                    x.store(&xa[part]);
                    xp.store(&xpa[part]);
                    y.store(&ya[part]);
                    yp.store(&ypa[part]);
                    cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_s_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);
                    FF_algorithm::thin_kicker_unit(xp, b_k);
                    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, reference_cdt * 0.5);

                    xa[part] = x;
                    xpa[part] = xp;
                    ya[part] = y;
                    ypa[part] = yp;
                    cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }

        }
        else
        {
            double b_k_pul = b_k / l;

            double step_reference_cdt = reference_cdt / steps;
            double step_length = length / steps;
            double step_strength[2] = { b_k_pul * step_length, 0.0 };

            // real particles
            {
                bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

                const int num_blocks = local_num / gsvsize;
                const int block_last = num_blocks * gsvsize;

                #pragma omp parallel for
                for (int part = 0; part < block_last; part += gsvsize)
                {
                    GSVector x(&xa[part]);
                    GSVector xp(&xpa[part]);
                    GSVector y(&ya[part]);
                    GSVector yp(&ypa[part]);
                    GSVector cdt(&cdta[part]);
                    GSVector dpop(&dpopa[part]);

                    FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_kicker_unit<GSVector>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                    x.store(&xa[part]);
                    xp.store(&xpa[part]);
                    y.store(&ya[part]);
                    yp.store(&ypa[part]);
                    cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::yoshida6<double, FF_algorithm::thin_kicker_unit<double>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                    xa[part] = x;
                    xpa[part] = xp;
                    ya[part] = y;
                    ypa[part] = yp;
                    cdta[part] = cdt;
                    dpopa[part] = dpop;
                }

            }

            // spectator particles
            {
                bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);

                const int num_blocks = local_s_num / gsvsize;
                const int block_last = num_blocks * gsvsize;

                #pragma omp parallel for
                for (int part = 0; part < block_last; part += gsvsize)
                {
                    GSVector x(&xa[part]);
                    GSVector xp(&xpa[part]);
                    GSVector y(&ya[part]);
                    GSVector yp(&ypa[part]);
                    GSVector cdt(&cdta[part]);
                    GSVector dpop(&dpopa[part]);

                    FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_kicker_unit<GSVector>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                    x.store(&xa[part]);
                    xp.store(&xpa[part]);
                    y.store(&ya[part]);
                    yp.store(&ypa[part]);
                    cdt.store(&cdta[part]);
                    dpop.store(&dpopa[part]);
                }

                for (int part = block_last; part < local_s_num; ++part)
                {
                    double x(xa[part]);
                    double xp(xpa[part]);
                    double y(ya[part]);
                    double yp(ypa[part]);
                    double cdt(cdta[part]);
                    double dpop(dpopa[part]);

                    FF_algorithm::yoshida6<double, FF_algorithm::thin_kicker_unit<double>, 1>
                        ( x, xp, y, yp, cdt, dpop,
                          pref, m, step_reference_cdt,
                          step_length, step_strength, steps );

                    xa[part] = x;
                    xpa[part] = xp;
                    ya[part] = y;
                    ypa[part] = yp;
                    cdta[part] = cdt;
                    dpopa[part] = dpop;
                }
            }
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
#endif
}

