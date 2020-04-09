#include "ff_sextupole.h"
#include "ff_algorithm.h"

namespace
{
    struct thin_kicker
    {
        Particles p;
        ParticleMasks m;

        double const* k;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { 
            if(m(i)) 
            {
                FF_algorithm::thin_sextupole_unit(
                        p(i, 0), p(i, 1), p(i, 2), p(i, 3), k); 
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
                FF_algorithm::yoshida<double, 
                    FF_algorithm::thin_sextupole_unit<double>, 2, 1>(
                            p(i, 0), p(i, 1), p(i, 2),
                            p(i, 3), p(i, 4), p(i, 5),
                            pref, mass, step_ref_cdt,
                            step_len, step_str, steps );
            }
        }
    };

    void apply_thin_kick(Bunch& bunch, ParticleGroup pg, 
            double const* k)
    {
        if(!bunch.get_local_num(pg)) return;

        auto parts = bunch.get_local_particles(pg);
        auto masks = bunch.get_local_particle_masks(pg);

        thin_kicker tk{parts, masks, k};
        Kokkos::parallel_for(bunch.size(pg), tk);
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

    void get_reference_cdt_zero(Reference_particle& ref, 
            double const* k)
    {
        // propagate the bunch design reference particle
        double x  = ref.get_state()[Bunch::x];
        double xp = ref.get_state()[Bunch::xp];
        double y  = ref.get_state()[Bunch::y];
        double yp = ref.get_state()[Bunch::yp];

        FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

        ref.set_state_xp(xp);
        ref.set_state_yp(yp);
        ref.set_state_cdt(0.0);
    }

    // non zero length
    double get_reference_cdt_yoshida(Reference_particle& ref, 
            double len, double const* k, int steps)
    {
        double pref = ref.get_momentum();
        double mass = ref.get_mass();
        double cdt = 0.0;

        auto  st = ref.get_state();

        // steps comes from base class, set in apply method
        double step_len = len / steps;

        // for >0 length, hk,vk is the strength/length of the kick
        double step_str[2] = { k[0]*step_len, k[1]*step_len };

        // propagate
        FF_algorithm::yoshida<double, 
            FF_algorithm::thin_sextupole_unit<double>, 4, 1>( 
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5],
                    pref, mass, 0.0,
                    step_len, step_str, steps );

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }

}


#if 0
double FF_sextupole::get_reference_cdt(double length, double * k, Reference_particle &reference_particle) 
{
    if (length == 0) 
    {
        reference_particle.set_state_cdt(0.0);
        return 0.0;
    } 
    else 
    {
        double reference_momentum = reference_particle.get_momentum();
        double m = reference_particle.get_mass();
        double step_length = length/steps;
        double step_strength[2] = { k[0]*step_length, k[1]*step_length };

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(0.0);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        FF_algorithm::yoshida<double, FF_algorithm::thin_sextupole_unit<double>, 4, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  0.0,
                  step_length, step_strength, steps );

        reference_particle.set_state(x, xp, y, yp, cdt, dpop);

        return cdt;
    }
}
#endif

void FF_sextupole::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error(
            "Propagate JetParticle through a sextupole element is yet to be implemented");
}

void FF_sextupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    // strength
    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k2");
    k[1] = slice.get_lattice_element().get_double_attribute("k2s", 0.0);

    // tilting
    double tilt = slice.get_lattice_element().get_double_attribute("tilt", 0.0);
    if (tilt != 0.0)
    {
        std::complex<double> ck2(k[0], +k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -3.0*tilt));
        k[0] = ck2.real();
        k[1] = ck2.imag();
    }

    // scaling
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
    double brho_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]) 
                                         / ref_l.get_charge();  // GV/c

    double scale = brho_l / brho_b;

    k[0] *= scale;
    k[1] *= scale;

    if ( close_to_zero(length) )
    {
        get_reference_cdt_zero(ref_l, k);

        apply_thin_kick(bunch, ParticleGroup::regular, k);
        apply_thin_kick(bunch, ParticleGroup::spectator, k);
    }
    else
    {
        double pref = bunch.get_reference_particle().get_momentum();
        double mass = bunch.get_mass();

        double ref_cdt = get_reference_cdt_yoshida(ref_l, length, k, steps);

        double step_ref_cdt = ref_cdt / steps;
        double step_len = length / steps;
        double step_str[2] = { k[0]*step_len, k[1]*step_len };

        // propagate
        apply_yoshida_kick(bunch, ParticleGroup::regular,
                pref, mass, step_ref_cdt, step_len, step_str, steps);

        apply_yoshida_kick(bunch, ParticleGroup::spectator,
                pref, mass, step_ref_cdt, step_len, step_str, steps);

        // trajectory
        bunch.get_reference_particle().increment_trajectory(length);
    }

#if 0
    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int local_num = bunch.get_local_num();
    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

    if (length == 0.0)
    {
        // propagate the bunch design reference particle
        double x  = ref_l.get_state()[Bunch::x];
        double xp = ref_l.get_state()[Bunch::xp];
        double y  = ref_l.get_state()[Bunch::y];
        double yp = ref_l.get_state()[Bunch::yp];
        double dpop = ref_l.get_state()[Bunch::dpop];

        FF_algorithm::thin_sextupole_unit(x, xp,  y, yp, k);
        ref_l.set_state(x, xp, y, yp, 0.0, dpop);

        // bunch particles
        {
            bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);
            const int local_num = bunch.get_local_num();

            const int gsvsize = GSVector::size();
            const int num_blocks = local_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector  x( &xa[part]);
                GSVector xp(&xpa[part]);
                GSVector  y( &ya[part]);
                GSVector yp(&ypa[part]);

                FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

                xp.store(&xpa[part]);
                yp.store(&ypa[part]);
            }

            for (int part = block_last; part < local_num; ++part)
            {
                double  x( xa[part]);
                double xp(xpa[part]);
                double  y( ya[part]);
                double yp(ypa[part]);

                FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

                xpa[part] = xp;
                ypa[part] = yp;
            }
        }

        // bunch spectator particles
        {
            bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);
            const int local_num = bunch.get_local_spectator_num();

            const int gsvsize = GSVector::size();
            const int num_blocks = local_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector  x( &xa[part]);
                GSVector xp(&xpa[part]);
                GSVector  y( &ya[part]);
                GSVector yp(&ypa[part]);

                FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

                xp.store(&xpa[part]);
                yp.store(&ypa[part]);
            }

            for (int part = block_last; part < local_num; ++part)
            {
                double  x( xa[part]);
                double xp(xpa[part]);
                double  y( ya[part]);
                double yp(ypa[part]);

                FF_algorithm::thin_sextupole_unit(x, xp, y, yp, k);

                xpa[part] = xp;
                ypa[part] = yp;
            }
        }
    }
    else
    {
        double reference_momentum = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k, ref_l);
        double step_reference_cdt = reference_cdt/steps;
        double step_length = length/steps;
        double step_strength[2] = { k[0]*step_length, k[1]*step_length };

        // bunch particles
        {
            bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);
            const int local_num = bunch.get_local_num();

            const int gsvsize = GSVector::size();
            const int num_blocks = local_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector    x(   &xa[part]);
                GSVector   xp(  &xpa[part]);
                GSVector    y(   &ya[part]);
                GSVector   yp(  &ypa[part]);
                GSVector  cdt( &cdta[part]);
                GSVector dpop(&dpopa[part]);

                FF_algorithm::yoshida<GSVector, FF_algorithm::thin_sextupole_unit<GSVector>, 2, 1 >
                        ( x, xp, y, yp, cdt, dpop,
                          reference_momentum, m,
                          step_reference_cdt,
                          step_length, step_strength, steps );

                   x.store(&xa[part]);
                  xp.store(&xpa[part]);
                   y.store(&ya[part]);
                  yp.store(&ypa[part]);
                 cdt.store(&cdta[part]);
            }

            for (int part = block_last; part < local_num; ++part)
            {
                double    x(   xa[part]);
                double   xp(  xpa[part]);
                double    y(   ya[part]);
                double   yp(  ypa[part]);
                double  cdt( cdta[part]);
                double dpop(dpopa[part]);

                FF_algorithm::yoshida<double, FF_algorithm::thin_sextupole_unit<double>, 2, 1 >
                        ( x, xp, y, yp, cdt, dpop,
                          reference_momentum, m,
                          step_reference_cdt,
                          step_length, step_strength, steps );

                   xa[part] = x;
                  xpa[part] = xp;
                   ya[part] = y;
                  ypa[part] = yp;
                 cdta[part] = cdt;
            }
        }

        // bunch spectator particles
        {
            bunch.set_spectator_arrays(xa, xpa, ya, ypa, cdta, dpopa);
            const int local_num = bunch.get_local_spectator_num();

            const int gsvsize = GSVector::size();
            const int num_blocks = local_num / gsvsize;
            const int block_last = num_blocks * gsvsize;

            #pragma omp parallel for
            for (int part = 0; part < block_last; part += gsvsize)
            {
                GSVector    x(   &xa[part]);
                GSVector   xp(  &xpa[part]);
                GSVector    y(   &ya[part]);
                GSVector   yp(  &ypa[part]);
                GSVector  cdt( &cdta[part]);
                GSVector dpop(&dpopa[part]);

                FF_algorithm::yoshida<GSVector, FF_algorithm::thin_sextupole_unit<GSVector>, 2, 1 >
                        ( x, xp, y, yp, cdt, dpop,
                          reference_momentum, m,
                          step_reference_cdt,
                          step_length, step_strength, steps );

                   x.store(&xa[part]);
                  xp.store(&xpa[part]);
                   y.store(&ya[part]);
                  yp.store(&ypa[part]);
                 cdt.store(&cdta[part]);
            }

            for (int part = block_last; part < local_num; ++part)
            {
                double    x(   xa[part]);
                double   xp(  xpa[part]);
                double    y(   ya[part]);
                double   yp(  ypa[part]);
                double  cdt( cdta[part]);
                double dpop(dpopa[part]);

                FF_algorithm::yoshida<double, FF_algorithm::thin_sextupole_unit<double>, 2, 1 >
                        ( x, xp, y, yp, cdt, dpop,
                          reference_momentum, m,
                          step_reference_cdt,
                          step_length, step_strength, steps );

                   xa[part] = x;
                  xpa[part] = xp;
                   ya[part] = y;
                  ypa[part] = yp;
                 cdta[part] = cdt;
            }
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
#endif
}

