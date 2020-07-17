#include "ff_elens.h"
#include "ff_algorithm.h"
#include "ff_patterned_propagator.h"
#include "synergia/utils/simple_timer.h"


void FF_elens::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error(
            "libFF::elens for JetParticles "
            "has yet to be implemented");
}

void FF_elens::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    scoped_simple_timer timer("libFF_elens");

    double length = slice.get_right() - slice.get_left();

    auto const& elem = slice.get_lattice_element();

    // hk and vk are the hk/vk under lattice reference momentum
    double   lengh = elem.get_double_attribute("l");
    double current = elem.get_double_attribute("current");
    double eenergy = elem.get_double_attribute("eenergy") * 0.001;
    double  radius = elem.get_double_attribute("radius");

    bool  gaussian = elem.get_double_attribute("gaussian", 0.0) != 0.0;
    bool   uniform = elem.get_double_attribute("uniform", 0.0) != 0.0;

    if (!(uniform || gaussian)) 
    {
        throw std::runtime_error(
                "elens must set either gaussian or uniform attribute");
    }

    if (gaussian && uniform) 
    {
        throw std::runtime_error(
                "elens must not set both gaussian and uniform attributes");
    }

    const double  current_over_e = current / pconstants::e;

    // references
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    const double gamma_e = (eenergy + pconstants::me) / pconstants::me;
    const double  beta_e = std::sqrt(1.0 - 1.0/(gamma_e*gamma_e));
    const double gamma_b = ref_b.get_gamma();
    const double  beta_b = ref_b.get_beta();

    const double  pref_l = ref_l.get_momentum();
    const double  pref_b = ref_b.get_momentum();

    const double     m_l = ref_l.get_mass();
    const double     m_b = ref_b.get_mass();

    const double k[6] = {
        beta_b, gamma_b, beta_e, current_over_e, length, radius
    };

    if (gaussian)
    {
        using pp = FF_patterned_propagator<double, 3,
              FF_algorithm::elens_kick_gaussian>;

        if (close_to_zero(length))
        {
            pp::get_reference_cdt_zero(ref_l, k);

            pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
            pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
        }
        else
        {
            double ref_cdt = pp::get_reference_cdt_simple(
                    ref_l, length, k);

            pp::apply_simple_kick(bunch, ParticleGroup::regular,
                    pref_b, m_b, ref_cdt, length, k);

            pp::apply_simple_kick(bunch, ParticleGroup::spectator,
                    pref_b, m_b, ref_cdt, length, k);
        }
    }
    else
    {
        using pp = FF_patterned_propagator<double, 3,
              FF_algorithm::elens_kick_uniform>;

        if (close_to_zero(length))
        {
            pp::get_reference_cdt_zero(ref_l, k);

            pp::apply_thin_kick(bunch, ParticleGroup::regular, k);
            pp::apply_thin_kick(bunch, ParticleGroup::spectator, k);
        }
        else
        {
            double ref_cdt = pp::get_reference_cdt_simple(
                    ref_l, length, k);

            pp::apply_simple_kick(bunch, ParticleGroup::regular,
                    pref_b, m_b, ref_cdt, length, k);

            pp::apply_simple_kick(bunch, ParticleGroup::spectator,
                    pref_b, m_b, ref_cdt, length, k);
        }

    }

    // increment the trajectory
    bunch.get_reference_particle().increment_trajectory(length);


#if 0
    double ref_x    = ref_l.get_state()[Bunch::x];
    double ref_xp   = ref_l.get_state()[Bunch::xp];
    double ref_y    = ref_l.get_state()[Bunch::y];
    double ref_yp   = ref_l.get_state()[Bunch::yp];
    double ref_cdt  = 0.0;
    double ref_dpop = ref_l.get_state()[Bunch::dpop];

    int local_num = bunch.get_local_num();
    int local_s_num = bunch.get_local_spectator_num();

    MArray2d_ref particles = bunch.get_local_particles();
    MArray2d_ref s_particles = bunch.get_local_spectator_particles();

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference particle which
        // relative to plattice.
        // also update the reference particle

        if (gaussian)
        {
            // reference particle
            FF_algorithm::elens_kick_gaussian(ref_x, ref_xp, ref_y, ref_yp, ref_dpop, k);
                    //beta_b, gamma_b, beta_e, current_over_e, length, radius);

            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::elens_kick_gaussian(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }

            // spectators
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x   (s_particles[part][Bunch::x   ]);
                double xp  (s_particles[part][Bunch::xp  ]);
                double y   (s_particles[part][Bunch::y   ]);
                double yp  (s_particles[part][Bunch::yp  ]);
                double cdt (s_particles[part][Bunch::cdt ]);
                double dpop(s_particles[part][Bunch::dpop]);

                FF_algorithm::elens_kick_gaussian(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                s_particles[part][Bunch::x]  = x;
                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::y]  = y;
                s_particles[part][Bunch::yp] = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }
        else
        {
            // reference particle
            FF_algorithm::elens_kick_uniform(ref_x, ref_xp, ref_y, ref_yp, ref_dpop, k);
                    //beta_b, gamma_b, beta_e, current_over_e, length, radius);

            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::elens_kick_uniform(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }

            // spectators
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x   (s_particles[part][Bunch::x   ]);
                double xp  (s_particles[part][Bunch::xp  ]);
                double y   (s_particles[part][Bunch::y   ]);
                double yp  (s_particles[part][Bunch::yp  ]);
                double cdt (s_particles[part][Bunch::cdt ]);
                double dpop(s_particles[part][Bunch::dpop]);

                FF_algorithm::elens_kick_uniform(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                s_particles[part][Bunch::x]  = x;
                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::y]  = y;
                s_particles[part][Bunch::yp] = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }
    }
    else
    {
        if (gaussian)
        {
            // reference particle
            FF_algorithm::drift_unit(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop, 
                    length * 0.5, pref_l, m_l, 0.0);

            FF_algorithm::elens_kick_gaussian(ref_x, ref_xp, ref_y, ref_yp, ref_dpop, k);
                    //beta_b, gamma_b, beta_e, current_over_e, length, radius);

            FF_algorithm::drift_unit(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop, 
                    length * 0.5, pref_l, m_l, 0.0);

            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                FF_algorithm::elens_kick_gaussian(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }

            // spectators
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x   (s_particles[part][Bunch::x   ]);
                double xp  (s_particles[part][Bunch::xp  ]);
                double y   (s_particles[part][Bunch::y   ]);
                double yp  (s_particles[part][Bunch::yp  ]);
                double cdt (s_particles[part][Bunch::cdt ]);
                double dpop(s_particles[part][Bunch::dpop]);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                FF_algorithm::elens_kick_gaussian(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                s_particles[part][Bunch::x]  = x;
                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::y]  = y;
                s_particles[part][Bunch::yp] = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }
        else
        {
            // reference particle
            FF_algorithm::drift_unit(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop, 
                    length * 0.5, pref_l, m_l, 0.0);

            FF_algorithm::elens_kick_uniform(ref_x, ref_xp, ref_y, ref_yp, ref_dpop, k);
                    //beta_b, gamma_b, beta_e, current_over_e, length, radius);

            FF_algorithm::drift_unit(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop, 
                    length * 0.5, pref_l, m_l, 0.0);

            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                FF_algorithm::elens_kick_uniform(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }

            // spectators
            #pragma omp parallel for
            for (int part = 0; part < local_s_num; ++part)
            {
                double x   (s_particles[part][Bunch::x   ]);
                double xp  (s_particles[part][Bunch::xp  ]);
                double y   (s_particles[part][Bunch::y   ]);
                double yp  (s_particles[part][Bunch::yp  ]);
                double cdt (s_particles[part][Bunch::cdt ]);
                double dpop(s_particles[part][Bunch::dpop]);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                FF_algorithm::elens_kick_uniform(x, xp, y, yp, dpop, k);
                        //beta_b, gamma_b, beta_e, current_over_e, length, radius);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                s_particles[part][Bunch::x]  = x;
                s_particles[part][Bunch::xp] = xp;
                s_particles[part][Bunch::y]  = y;
                s_particles[part][Bunch::yp] = yp;
                s_particles[part][Bunch::cdt] = cdt;
            }
        }
    }

    // update reference particles state
    ref_l.set_state(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop);

    // increment the trajectory
    bunch.get_reference_particle().increment_trajectory(length);
#endif
}

