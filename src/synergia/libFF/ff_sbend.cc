#include "ff_sbend.h"
#include "synergia/lattice/chef_utils.h"

FF_sbend::FF_sbend()
{

}

double FF_sbend::get_reference_cdt(double length, double strength, double angle, 
                                   bool ledge, bool redge,
                                   double e1, double e2, double dphi,
                                   std::complex<double> const & phase,
                                   std::complex<double> const & term,
                                   Reference_particle &reference_particle) 
{
    double reference_cdt;

    if (length == 0) 
    {
        reference_cdt = 0.0;
    } 
    else 
    {
        double pref = reference_particle.get_momentum();
        double m = reference_particle.get_mass();

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(reference_particle.get_state()[Bunch::cdt]);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double ce1 = cos(-e1);
        double se1 = sin(-e1);
        double ce2 = cos(-e2);
        double se2 = sin(-e2);

        double cdt_orig = cdt;

        if (ledge)
        {
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref, m);
            //FF_algorithm::edge_unit(y, yp, us_edge_k);
         }

        FF_algorithm::bend_complete(x, xp, y, yp, cdt, dpop,
                   dphi, strength, pref, m, 0.0/*ref cdt*/, phase, term);

        if (redge)
        {
            //FF_algorithm::edge_unit(y, yp, ds_edge_k);
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref, m);
        }

        reference_cdt = cdt - cdt_orig;
    }

    return reference_cdt;
}


double FF_sbend::get_reference_cdt(double length, double angle, double strength,
                                   Reference_particle &reference_particle) 
{
    double reference_cdt;
    if (length == 0) {
        reference_cdt = 0.0;
    } else {
        double reference_momentum = reference_particle.get_momentum();
        double m = reference_particle.get_mass();

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(reference_particle.get_state()[Bunch::cdt]);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        std::complex<double> phase = std::exp( std::complex<double>(0.0, angle) );
        std::complex<double> term = std::complex<double>(0.0, length / angle) *
                                    std::complex<double>(1.0 - cos(angle), - sin(angle));

        double cdt_orig = cdt;
        sbend_unit2(x, xp, y, yp, cdt, dpop,
                    length, angle, strength,
                    reference_momentum, m,
                    0.0, phase, term );
        reference_cdt = cdt - cdt_orig;
    }
    return reference_cdt;
}

void FF_sbend::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF sbend on JetParticle not implemented");

#if 0
    double length = slice.get_right() - slice.get_left();
    double angle = slice.get_lattice_element().get_double_attribute("angle");

    double cos_angle = cos(angle);
    double sin_angle = sin(angle);

    typedef PropagatorTraits<JetParticle>::State_t State_t;
    typedef PropagatorTraits<JetParticle>::Component_t Component_t;

    State_t& state = jet_particle.State();

    Component_t & x(state[Chef::x]);
    Component_t & xp(state[Chef::xp]);
    Component_t & y(state[Chef::y]);
    Component_t & yp(state[Chef::yp]);
    Component_t & cdt(state[Chef::cdt]);
    Component_t & dpop(state[Chef::dpop]);

    double reference_momentum = jet_particle.ReferenceMomentum();
    double reference_brho     = jet_particle.ReferenceBRho();
    double m = jet_particle.Mass();

    sbend_unit(x, xp, y, yp, cdt, dpop,
               length, cos_angle, sin_angle,
               reference_momentum, m, reference_brho);
#endif
}

void FF_sbend::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double      a = slice.get_lattice_element().get_double_attribute("angle");
    double      l = slice.get_lattice_element().get_double_attribute("l");
    double length = slice.get_right() - slice.get_left();
    double  angle = ( length / l ) * a;

    double     r0 = l / a;

    double ledge  = slice.has_left_edge();
    double redge  = slice.has_right_edge();

    double cos_angle = cos(angle);
    double sin_angle = sin(angle);

    double e1 = 0.0;
    double e2 = 0.0;

    if (slice.get_lattice_element().has_double_attribute("e1"))
        e1 = slice.get_lattice_element().get_double_attribute("e1");
    if (slice.get_lattice_element().has_double_attribute("e2"))
        e2 = slice.get_lattice_element().get_double_attribute("e2");


    int cf = 0;  // combined function
    double kl[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    if (slice.get_lattice_element().has_double_attribute("k1"))
    {
        // quad component
        kl[0] = slice.get_lattice_element().get_double_attribute("k1");
        kl[1] = 0.0;

        if (kl[0] != 0) cf = 1;
    }

    if (slice.get_lattice_element().has_double_attribute("k2"))
    {
        // sextupole component
        kl[2] = slice.get_lattice_element().get_double_attribute("k2");
        kl[3] = 0.0;

        if (kl[2] != 0) cf = 2;
    }

    if (slice.get_lattice_element().has_double_attribute("k3"))
    {
        // octupole component
        kl[4] = slice.get_lattice_element().get_double_attribute("k3");
        kl[5] = 0.0;

        if (kl[4] != 0) cf = 3;
    }


    double usAngle = e1;
    double dsAngle = -e2;
    double usFaceAngle = e1;
    double dsFaceAngle = e2;

    if (!redge)
    {
        dsAngle = 0;
        dsFaceAngle = 0;
    }

    if (!ledge)
    {
        usAngle = 0;
        usFaceAngle = 0;
    }

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double reference_momentum = bunch.get_reference_particle().get_momentum();
    double reference_brho     = reference_momentum / PH_CNV_brho_to_p;
    int    reference_charge   = bunch.get_reference_particle().get_charge();
    double m = bunch.get_mass();

    double strength = reference_brho * a / l;

    double psi = angle - (usFaceAngle + dsFaceAngle);
    double dphi = -psi;
    std::complex<double> phase = std::exp( std::complex<double>(0.0, psi) );
    std::complex<double> term = std::complex<double>(0.0, length / angle) *
                                std::complex<double>(1.0 - cos_angle, - sin_angle) *
                                std::complex<double>(cos(dsFaceAngle), -sin(dsFaceAngle));

    double pref = reference_momentum;

    double ce1 = cos(-e1);
    double se1 = sin(-e1);
    double ce2 = cos(-e2);
    double se2 = sin(-e2);

    double us_edge_k =   ((reference_charge > 0) ? 1.0 : -1.0) * strength * tan(usAngle) / reference_brho;
    double ds_edge_k = - ((reference_charge > 0) ? 1.0 : -1.0) * strength * tan(dsAngle) / reference_brho;

    double us_edge_k_p =   ((reference_charge > 0) ? 1.0 : -1.0) * strength / reference_brho;
    double ds_edge_k_p = - ((reference_charge > 0) ? 1.0 : -1.0) * strength / reference_brho;

    double reference_cdt = get_reference_cdt(length, strength, angle, ledge, redge, e1, e2, dphi,
            phase, term, bunch.get_reference_particle());

    std::complex<double> phase_e1 = FF_algorithm::bend_edge_phase(e1);
    std::complex<double> phase_e2 = FF_algorithm::bend_edge_phase(e2);

    std::complex<double> phase_unit = FF_algorithm::bend_unit_phase(angle);
    std::complex<double>  term_unit = FF_algorithm::bend_unit_term(r0, angle);

    double step_reference_cdt = reference_cdt / steps;
    double step_angle = angle / steps;
    double step_length = length / steps;
    double step_strength[6] = { kl[0] * step_length, kl[1] * step_length,
                                kl[2] * step_length, kl[3] * step_length,
                                kl[2] * step_length, kl[3] * step_length };

    if (cf == 0)
    {
        // no combined function
        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part) 
        {
            double x   (particles[part][Bunch::x   ]);
            double xp  (particles[part][Bunch::xp  ]);
            double y   (particles[part][Bunch::y   ]);
            double yp  (particles[part][Bunch::yp  ]);
            double cdt (particles[part][Bunch::cdt ]);
            double dpop(particles[part][Bunch::dpop]);

            if (ledge)
            {
                // slot
                FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref, m);

                // edge
                // FF_algorithm::edge_unit(y, yp, us_edge_k);
                FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);
            }

            // bend
            FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop,
                       dphi, strength, pref, m, reference_cdt,
                       phase, term);

            if (redge)
            {
                // edge
                // FF_algorithm::edge_unit(y, yp, ds_edge_k);
                FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);

                // slot
                FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref, m);
            }

            particles[part][Bunch::x]  = x;
            particles[part][Bunch::xp] = xp;
            particles[part][Bunch::y]  = y;
            particles[part][Bunch::yp] = yp;
            particles[part][Bunch::cdt] = cdt;
        }

    }
    else
    {
        // with combined function
        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part) 
        {
            double x   (particles[part][Bunch::x   ]);
            double xp  (particles[part][Bunch::xp  ]);
            double y   (particles[part][Bunch::y   ]);
            double yp  (particles[part][Bunch::yp  ]);
            double cdt (particles[part][Bunch::cdt ]);
            double dpop(particles[part][Bunch::dpop]);

            if (ledge)
            {
                FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce1, se1, pref, m);

                //FF_algorithm::edge_unit(y, yp, us_edge_k);
                FF_algorithm::edge_unit(y, xp, yp, dpop, us_edge_k_p);

                FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e1, phase_e1, strength, pref, m);
            }

            // bend
            FF_algorithm::bend_yoshida4<double, FF_algorithm::thin_cf_kick_2<double>, 2>
                ( x, xp, y, yp, cdt, dpop,
                  pref, m, step_reference_cdt,
                  step_angle, step_strength,
                  r0, strength, steps );

            if (redge)
            {
                FF_algorithm::bend_edge(x, xp, y, yp, cdt, dpop, e2, phase_e2, strength, pref, m);

                //FF_algorithm::edge_unit(y, yp, ds_edge_k);
                FF_algorithm::edge_unit(y, xp, yp, dpop, ds_edge_k_p);

                FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ce2, se2, pref, m);
            }

            particles[part][Bunch::x]  = x;
            particles[part][Bunch::xp] = xp;
            particles[part][Bunch::y]  = y;
            particles[part][Bunch::yp] = yp;
            particles[part][Bunch::cdt] = cdt;
        }
    }

    bunch.get_reference_particle().increment_trajectory(length);
}

template<class Archive>
    void
    FF_sbend::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_sbend::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_sbend::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_sbend::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_sbend::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_sbend::~FF_sbend()
{

}

