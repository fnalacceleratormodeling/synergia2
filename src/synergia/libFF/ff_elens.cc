#include "ff_elens.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

double FF_elens::get_reference_cdt(double length, double k, Reference_particle &reference_particle, bool simple)
{
    double pref = reference_particle.get_momentum();
    double m = reference_particle.get_mass();

    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    if ( close_to_zero(length) )
    {
        // for 0 length, hk,vk is the total strength of the kick
        FF_algorithm::thin_kicker_unit(xp, k);
    }
    else
    {
        if (simple)
        {
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, 0.0);

            // for >0 length, hk,vk is the strength/length of the kick
            FF_algorithm::thin_kicker_unit(xp, k * length);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length * 0.5, pref, m, 0.0);
        }
        else
        {
            // steps comes from base class, set in apply method
            double step_length = length / steps;

            // for >0 length, hk,vk is the strength/length of the kick
            double step_strength[2] = { k*step_length, 0.0 };

            // propagate
            FF_algorithm::yoshida6<double, FF_algorithm::thin_kicker_unit<double>, 1>
                ( x, xp, y, yp, cdt, dpop,
                  pref, m, 0.0,
                  step_length, step_strength, steps );
        }
    }

    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}

void FF_elens::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF::elens for JetParticles has yet to be implemented");
}

void FF_elens::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    // hk and vk are the hk/vk under lattice reference momentum
    double   lengh = slice.get_lattice_element().get_double_attribute("l");
    double current = slice.get_lattice_element().get_double_attribute("current");
    double eenergy = slice.get_lattice_element().get_double_attribute("eenergy") * 0.001;
    double  radius = slice.get_lattice_element().get_double_attribute("radius");

    bool  gaussian = slice.get_lattice_element().get_double_attribute("gaussian") != 0.0;
    bool   uniform = slice.get_lattice_element().get_double_attribute("uniform") != 0.0;

    if (!(uniform || gaussian)) 
    {
        throw std::runtime_error("elens must set either gaussian or uniform attribute");
    }

    if (gaussian && uniform) 
    {
        throw std::runtime_error("elens must not set both gaussian and uniform attributes");
    }

    const double  current_over_e = current / pconstants::e;
    const double    elens_length = length;
    const double electron_energy = eenergy;

    // references
    Reference_particle       & ref_l = bunch.get_design_reference_particle();
    Reference_particle const & ref_b = bunch.get_reference_particle();

    double ref_x    = ref_l.get_state()[Bunch::x];
    double ref_xp   = ref_l.get_state()[Bunch::xp];
    double ref_y    = ref_l.get_state()[Bunch::y];
    double ref_yp   = ref_l.get_state()[Bunch::yp];
    double ref_cdt  = 0.0;
    double ref_dpop = ref_l.get_state()[Bunch::dpop];

    const double gamma_e = (eenergy + pconstants::me) / pconstants::me;
    const double  beta_e = std::sqrt(1.0 - 1.0/(gamma_e*gamma_e));
    const double gamma_b = ref_b.get_gamma();
    const double  beta_b = ref_b.get_beta();

    const double  pref_l = ref_l.get_momentum();
    const double  pref_b = ref_b.get_momentum();

    const double     m_l = ref_l.get_mass();
    const double     m_b = ref_b.get_mass();

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference particle which
        // relative to plattice.
        // also update the reference particle

        if (gaussian)
        {
            // reference particle
            FF_algorithm::elens_kick_gaussian(ref_x, ref_xp, ref_y, ref_yp, ref_dpop,
                    beta_b, gamma_b, beta_e, current_over_e, length, radius);

            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::elens_kick_gaussian(x, xp, y, yp, dpop,
                        beta_b, gamma_b, beta_e, current_over_e, length, radius);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }
        }
        else
        {
            // reference particle
            FF_algorithm::elens_kick_uniform(ref_x, ref_xp, ref_y, ref_yp, ref_dpop,
                    beta_b, gamma_b, beta_e, current_over_e, length, radius);

            #pragma omp parallel for
            for (int part = 0; part < local_num; ++part)
            {
                double x   (particles[part][Bunch::x   ]);
                double xp  (particles[part][Bunch::xp  ]);
                double y   (particles[part][Bunch::y   ]);
                double yp  (particles[part][Bunch::yp  ]);
                double cdt (particles[part][Bunch::cdt ]);
                double dpop(particles[part][Bunch::dpop]);

                FF_algorithm::elens_kick_uniform(x, xp, y, yp, dpop,
                        beta_b, gamma_b, beta_e, current_over_e, length, radius);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
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

            FF_algorithm::elens_kick_gaussian(ref_x, ref_xp, ref_y, ref_yp, ref_dpop,
                    beta_b, gamma_b, beta_e, current_over_e, length, radius);

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

                FF_algorithm::elens_kick_gaussian(x, xp, y, yp, dpop,
                        beta_b, gamma_b, beta_e, current_over_e, length, radius);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }
        }
        else
        {
            // reference particle
            FF_algorithm::drift_unit(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop, 
                    length * 0.5, pref_l, m_l, 0.0);

            FF_algorithm::elens_kick_uniform(ref_x, ref_xp, ref_y, ref_yp, ref_dpop,
                    beta_b, gamma_b, beta_e, current_over_e, length, radius);

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

                FF_algorithm::elens_kick_uniform(x, xp, y, yp, dpop,
                        beta_b, gamma_b, beta_e, current_over_e, length, radius);

                FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, 
                        length * 0.5, pref_b, m_b, ref_cdt * 0.5);

                particles[part][Bunch::x]  = x;
                particles[part][Bunch::xp] = xp;
                particles[part][Bunch::y]  = y;
                particles[part][Bunch::yp] = yp;
                particles[part][Bunch::cdt] = cdt;
            }
        }
    }

    // update reference particles state
    ref_l.set_state(ref_x, ref_xp, ref_y, ref_yp, ref_cdt, ref_dpop);

    // increment the trajectory
    bunch.get_reference_particle().increment_trajectory(length);
}

template<class Archive>
    void
    FF_elens::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_elens::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_elens::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_elens::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_elens::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_elens::~FF_elens()
{

}

