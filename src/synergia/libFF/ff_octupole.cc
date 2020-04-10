#include "ff_octupole.h"
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
                FF_algorithm::thin_octupole_unit(
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
                    FF_algorithm::thin_octupole_unit<double>, 4, 1>(
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

        FF_algorithm::thin_octupole_unit(x, xp, y, yp, k);

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
            FF_algorithm::thin_octupole_unit<double>, 4, 1>( 
                    st[0], st[1], st[2], 
                    st[3], cdt, st[5],
                    pref, mass, 0.0,
                    step_len, step_str, steps );

        st[4] = cdt;
        ref.set_state(st);
        return cdt;
    }

}


void FF_octupole::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error(
            "Propagate JetParticle through a octupole element is yet to be implemented");
}

void FF_octupole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    double k[2];
    k[0] = slice.get_lattice_element().get_double_attribute("k3", 0.0);
    k[1] = slice.get_lattice_element().get_double_attribute("k3s", 0.0);

    // tilting
    double tilt = slice.get_lattice_element().get_double_attribute("tilt");
    if (tilt != 0.0)
    {
        std::complex<double> ck2(k[0], -k[1]);
        ck2 = ck2 * exp(std::complex<double>(0.0, -4.0*tilt));
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
    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    if (length == 0.0)
    {
        // propagate the bunch design reference particle
        double x  = ref_l.get_state()[Bunch::x];
        double xp = ref_l.get_state()[Bunch::xp];
        double y  = ref_l.get_state()[Bunch::y];
        double yp = ref_l.get_state()[Bunch::yp];
        double dpop = ref_l.get_state()[Bunch::dpop];

        FF_algorithm::thin_octupole_unit(x, xp,  y, yp, k);
        ref_l.set_state(x, xp, y, yp, 0.0, dpop);

        // bunch particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x(particles[part][Bunch::x]);
                double xp(particles[part][Bunch::xp]);
                double y(particles[part][Bunch::y]);
                double yp(particles[part][Bunch::yp]);

                FF_algorithm::thin_octupole_unit(x, xp, y, yp, k);

                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::yp] = yp;
           }
        }

        // bunch spectator particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x(s_particles[part][Bunch::x]);
                double xp(s_particles[part][Bunch::xp]);
                double y(s_particles[part][Bunch::y]);
                double yp(s_particles[part][Bunch::yp]);

                FF_algorithm::thin_octupole_unit(x, xp, y, yp, k);

                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::yp] = yp;
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
            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x(particles[part][Bunch::x]);
                double xp(particles[part][Bunch::xp]);
                double y(particles[part][Bunch::y]);
                double yp(particles[part][Bunch::yp]);
                double cdt(particles[part][Bunch::cdt]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::yoshida<double, FF_algorithm::thin_octupole_unit<double>, 4, 1 >
                        ( x, xp, y, yp, cdt, dpop,
                          reference_momentum, m,
                          step_reference_cdt,
                          step_length, step_strength, steps );

                particles[part][Bunch::x] = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y] = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }
        }

        // bunch spectator particles
        {
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x(s_particles[part][Bunch::x]);
                double xp(s_particles[part][Bunch::xp]);
                double y(s_particles[part][Bunch::y]);
                double yp(s_particles[part][Bunch::yp]);
                double cdt(s_particles[part][Bunch::cdt]);
                double dpop(s_particles[part][Bunch::dpop]);

                FF_algorithm::yoshida<double, FF_algorithm::thin_octupole_unit<double>, 4, 1 >
                        ( x, xp, y, yp, cdt, dpop,
                          reference_momentum, m,
                          step_reference_cdt,
                          step_length, step_strength, steps );

                s_particles[part][Bunch::x] = x;
                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::y] = y;
                s_particles[part][Bunch::yp] = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
#endif
}

template<class Archive>
    void
    FF_octupole::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_octupole::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_octupole::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_octupole::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_octupole::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_octupole::~FF_octupole()
{

}

