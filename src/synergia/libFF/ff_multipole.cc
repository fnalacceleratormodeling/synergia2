
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
        ConstParticleMasks masks;
        const MultipoleParams mp;

        KOKKOS_INLINE_FUNCTION
        void operator()(const int i) const
        {
            if (masks(i))
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
                                p(i,0), p(i,1), p(i,2), p(i,3), &mp.kl[n*2], n+1);
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

    auto const& element = slice.get_lattice_element();

    // extract attributes
    if ( element.has_vector_attribute("knl") 
            || element.has_vector_attribute("ksl") )
    {
        // it is in Mad X format
        std::vector<double> k0(1, 0.0);

        knl = element.get_vector_attribute("knl", k0);
        ksl = element.get_vector_attribute("ksl", k0);

        if (knl.size() > ksl.size()) ksl.resize(knl.size(), 0.0);
        else if (knl.size() < ksl.size()) knl.resize(ksl.size(), 0.0);

        double tilt = element.get_double_attribute("tilt", 0.0);
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

            knl.push_back( element.get_double_attribute(skn, 0.0) );
             tn.push_back( element.get_double_attribute(stn, 0.0) );
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
    int num = bunch.size(ParticleGroup::regular);
    auto parts = bunch.get_local_particles(ParticleGroup::regular);
    auto masks = bunch.get_local_particle_masks(ParticleGroup::regular);

    PropMultipole multipole{parts, masks, mp};
    Kokkos::parallel_for(num, multipole);
}

