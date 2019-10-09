
#include "ff_multipole.h"
#include "ff_algorithm.h"

#include "synergia/utils/simple_timer.h"


namespace
{
    const int max_order = 8;

    struct MultipoleParams
    {
        bool kn[max_order];
        double kl[max_order*2];
    };

    void zero_params(MultipoleParams & mp)
    {
        for(int i=0; i<max_order; ++i)
        {
            mp.kn[i] = false;
            mp.kl[i*2+0] = 0.0;
            mp.kl[i*2+1] = 0.0;
        }
    }

    struct PropMultipole
    {
        Particles p;
        const MultipoleParams mp;
        const_k1b_dev valid;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (valid(i))
            {
                if (mp.kn[0])
                FF_algorithm::thin_dipole_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), &mp.kl[0]);

                if (mp.kn[1])
                FF_algorithm::thin_quadrupole_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), &mp.kl[2]);

                if (mp.kn[2])
                FF_algorithm::thin_sextupole_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), &mp.kl[4]);

                if (mp.kn[3])
                FF_algorithm::thin_octupole_unit(
                        p(i,0), p(i,1), p(i,2), p(i,3), &mp.kl[6]);

#if 1
                for(int n=4; n<max_order; ++n)
                {
                    if (mp.kn[n])
                        FF_algorithm::thin_magnet_unit(
                                p(i,0), p(i,1), p(i,2), p(i,3),
                                &mp.kl[n*2], n+1);
                }
#endif
            }
        }
    };
}


void FF_multipole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    scoped_simple_timer timer("libFF_multipole");

    if (slice.get_right() - slice.get_left() > 0.0)
        throw std::runtime_error("FF_multipole::apply() cannot deal with thick elements");

    MultipoleParams mp;
    zero_params(mp);

    std::vector<double> knl;
    std::vector<double> ksl;
    std::vector<double> tn;

    // extract attributes
    if ( slice.get_lattice_element().has_vector_attribute("knl") )
    {
        // it is in Mad X format
        std::vector<double> k0(1, 0.0);

        knl = slice.get_lattice_element().get_vector_attribute("knl");
        ksl = slice.get_lattice_element().get_vector_attribute("ksl", k0);

        if (knl.size() > ksl.size()) ksl.resize(knl.size(), 0.0);
        else if (knl.size() < ksl.size()) knl.resize(ksl.size(), 0.0);

        double tilt = slice.get_lattice_element().get_double_attribute("tilt", 0.0);
        tn.resize(knl.size(), tilt);
    }
    else
    {
        // in Mad 8 format
        std::string skn("k0l");
        std::string stn("t0");

        for (int i=0; i<6; ++i)
        {
            skn[1] = '0' + i;
            stn[1] = '0' + i;

            knl.push_back( slice.get_lattice_element().get_double_attribute(skn, 0.0) );
             tn.push_back( slice.get_lattice_element().get_double_attribute(stn, 0.0) );
        }

        int tail = knl.size()-1;
        while (tail && knl[tail] == 0.0) --tail;

        knl.resize(tail+1);
        ksl.resize(knl.size(), 0.0);
         tn.resize(knl.size(), 0.0);
    }

    // tilting
    for (int i=0; i<knl.size(); ++i)
    {
        if (tn[i] != 0.0)
        {
            std::complex<double> ck2(knl[i], ksl[i]);
            ck2 = ck2 * exp(std::complex<double>(0.0, -(i+1)*tn[i]));
            knl[i] = ck2.real();
            ksl[i] = ck2.imag();
        }
    }

    // scaling
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    double brho_l = ref_l.get_momentum() / ref_l.get_charge();  // GV/c
    double brho_b = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]) / ref_l.get_charge();  // GV/c

    // scale
    double scale = brho_l / brho_b;
    for (int i=0; i<knl.size(); ++i)
    {
        knl[i] *= scale;
        ksl[i] *= scale;
    }

    // prepare the params
    for (int i=0; i<knl.size(); ++i)
    {
        if (knl[i] || ksl[i]) mp.kn[i] = true;

        mp.kl[i*2+0] = knl[i];
        mp.kl[i*2+1] = ksl[i];
    }

    // propagate and update the design reference particle
    double x  = ref_l.get_state()[Bunch::x];
    double xp = ref_l.get_state()[Bunch::xp];
    double y  = ref_l.get_state()[Bunch::y];
    double yp = ref_l.get_state()[Bunch::yp];
    double dpop = ref_l.get_state()[Bunch::dpop];

    for (int n = 0; n < max_order; ++n) {
        if (mp.kn[n]) {
            FF_algorithm::thin_magnet_unit(x, xp, y, yp, &mp.kl[n*2], n+1);
        }
    }

    ref_l.set_state(x, xp, y, yp, 0.0, dpop);

    // bunch particles
    int num = bunch.get_local_num_slots(ParticleGroup::regular);
    auto parts = bunch.get_local_particles(ParticleGroup::regular);
    auto valid = bunch.get_local_particles_valid(ParticleGroup::regular);

    PropMultipole multipole{parts, mp, valid};
    Kokkos::parallel_for(num, multipole);

#if 0
    // propagate bunch particles
    for (int part = 0; part < local_num; ++part) 
    {
        double x   (particles[part][Bunch::x   ]);
        double xp  (particles[part][Bunch::xp  ]);
        double y   (particles[part][Bunch::y   ]);
        double yp  (particles[part][Bunch::yp  ]);

#if 1
        // dipole
        if (knl.size() > 0 && (knl[0] || ksl[0])) 
        {
            kL[0] = knl[0]; kL[1] = ksl[0];
            FF_algorithm::thin_dipole_unit(x, xp, y, yp, kL);
        }

        // quad
        if (knl.size() > 1 && (knl[1] || ksl[1])) 
        {
            kL[0] = knl[1]; kL[1] = ksl[1];
            FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kL);
        }

        // sextu
        if (knl.size() > 2 && (knl[2] || ksl[2])) 
        {
            kL[0] = knl[2]; kL[1] = ksl[2];
            FF_algorithm::thin_sextupole_unit(x, xp, y, yp, kL);
        }

        // octu 
        if (knl.size() > 3 && (knl[3] || ksl[3])) 
        {
            kL[0] = knl[3]; kL[1] = ksl[3];
            FF_algorithm::thin_octupole_unit(x, xp, y, yp, kL);
        }
#endif

        // higher orders
        for (int n = 4; n < knl.size(); ++n) {
            if (knl[n] || ksl[n]) {
                kL[0] = knl[n]; kL[1] = ksl[n];
                FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL, n+1);
            }
        }

        particles[part][Bunch::xp] = xp;
        particles[part][Bunch::yp] = yp;
    }

    // propagate bunch spectator particles
    for (int part = 0; part < local_s_num; ++part) 
    {
        double x   (s_particles[part][Bunch::x   ]);
        double xp  (s_particles[part][Bunch::xp  ]);
        double y   (s_particles[part][Bunch::y   ]);
        double yp  (s_particles[part][Bunch::yp  ]);

#if 1
        // dipole
        if (knl.size() > 0 && (knl[0] || ksl[0])) 
        {
            kL[0] = knl[0]; kL[1] = ksl[0];
            FF_algorithm::thin_dipole_unit(x, xp, y, yp, kL);
        }

        // quad
        if (knl.size() > 1 && (knl[1] || ksl[1])) 
        {
            kL[0] = knl[1]; kL[1] = ksl[1];
            FF_algorithm::thin_quadrupole_unit(x, xp, y, yp, kL);
        }

        // sextu
        if (knl.size() > 2 && (knl[2] || ksl[2])) 
        {
            kL[0] = knl[2]; kL[1] = ksl[2];
            FF_algorithm::thin_sextupole_unit(x, xp, y, yp, kL);
        }

        // octu 
        if (knl.size() > 3 && (knl[3] || ksl[3])) 
        {
            kL[0] = knl[3]; kL[1] = ksl[3];
            FF_algorithm::thin_octupole_unit(x, xp, y, yp, kL);
        }
#endif

        // higher orders
        for (int n = 4; n < knl.size(); ++n) {
            if (knl[n] || ksl[n]) {
                kL[0] = knl[n]; kL[1] = ksl[n];
                FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL, n+1);
            }
        }

        s_particles[part][Bunch::xp] = xp;
        s_particles[part][Bunch::yp] = yp;
    }
#endif
}

