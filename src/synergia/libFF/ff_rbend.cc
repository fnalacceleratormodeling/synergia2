#include "ff_rbend.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

FF_rbend::FF_rbend()
{

}

double FF_rbend::get_reference_cdt(double length, double strength, double angle,
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

        double theta = angle / 2.0;
        double ct = cos(-theta);
        double st = sin(-theta);

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(reference_particle.get_state()[Bunch::cdt]);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double cdt_orig = cdt;

        FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

        FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop, 
                    0.0, strength, pref, m, 0.0, phase, term);

        FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

        reference_cdt = cdt - cdt_orig;
    }

    return reference_cdt;
}


double FF_rbend::get_reference_cdt(double length, double * k,
                                   Reference_particle &reference_particle) 
{
    double reference_cdt;
    if (length == 0) {
        reference_cdt = 0.0;
    } else {
        double reference_momentum = reference_particle.get_momentum();
        double m = reference_particle.get_mass();
        double step_length = length/steps;
        double step_strength[6] = 
            { k[0]*step_length, k[1]*step_length,
              k[2]*step_length, k[3]*step_length,
              k[4]*step_length, k[5]*step_length };

        double x(reference_particle.get_state()[Bunch::x]);
        double xp(reference_particle.get_state()[Bunch::xp]);
        double y(reference_particle.get_state()[Bunch::y]);
        double yp(reference_particle.get_state()[Bunch::yp]);
        double cdt(reference_particle.get_state()[Bunch::cdt]);
        double dpop(reference_particle.get_state()[Bunch::dpop]);

        double cdt_orig = cdt;

        FF_algorithm::yoshida<double, FF_algorithm::thin_rbend_unit<double>, 4, 3 >
            ( x, xp, y, yp, cdt, dpop,
              reference_momentum, m,
              0.0, step_length, step_strength, steps );

        reference_cdt = cdt - cdt_orig;
    }
    return reference_cdt;
}

void FF_rbend::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF rbend on JetParticle not implemented");

#if 0
    double length = slice.get_right() - slice.get_left();
    double angle = slice.get_lattice_element().get_double_attribute("angle");
    double l     = slice.get_lattice_element().get_double_attribute("l");

    double k[6];
    k[2]  = slice.get_lattice_element().get_double_attribute("k1");
    k[3]  = 0;  // k1s
    k[4]  = slice.get_lattice_element().get_double_attribute("k2");
    k[5]  = 0;  // k2s

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

    k[0] = reference_brho * angle / l;
    k[1] = 0;
#endif
}

void FF_rbend::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double arc_length = slice.get_right() - slice.get_left();

    double  angle = slice.get_lattice_element().get_double_attribute("angle");
    double      l = slice.get_lattice_element().get_double_attribute("l");
    double    rho = l / ( 2.0 * sin( angle / 2.0 ) );
    double    arc = angle * rho;
    double length = l * (arc_length / arc);

    double k[6];

    k[0] = 2.0 * sin( angle / 2.0 ) / l;
    k[1] = 0;

    k[2]  = slice.get_lattice_element().get_double_attribute("k1");
    k[3]  = 0;  // k1s

    k[4]  = slice.get_lattice_element().get_double_attribute("k2");
    k[5]  = 0;  // k2s

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double reference_momentum = bunch.get_reference_particle().get_momentum();
    int    reference_charge   = bunch.get_reference_particle().get_charge();
    double m = bunch.get_mass();

    double pref = reference_momentum;
    double theta = angle / 2.0;
    double ct = cos(-theta);
    double st = sin(-theta);

    double reference_brho = pref / PH_CNV_brho_to_p;
    double strength = reference_brho * angle / arc;

    double usFaceAngle = 0.0;  // always 0 for rbends
    double dsFaceAngle = 0.0;  // alwasy 0 for rbends
    double psi = - ( usFaceAngle + dsFaceAngle );
    double dphi = -psi;
    std::complex<double> phase = std::exp( std::complex<double>(0.0, psi) );
    std::complex<double> term = length * std::complex<double> ( cos(dsFaceAngle), -sin(dsFaceAngle) );

    double edge_k = ((reference_charge > 0) ? 1.0 : -1.0 ) * strength * tan(theta) / reference_brho;

    if (k[2] == 0.0 && k[4] == 0.0)
    {
        double reference_cdt = get_reference_cdt(length, strength, angle, phase, term, 
                bunch.get_reference_particle());

        double * RESTRICT xa, * RESTRICT xpa;
        double * RESTRICT ya, * RESTRICT ypa;
        double * RESTRICT cdta, * RESTRICT dpopa;

        bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

        const int num_blocks = local_num / GSVector::size;
        const int block_last = num_blocks * GSVector::size;

#if 0
        // use the exact solution for dipole
        for (int part = 0; part < block_last; part += GSVector::size)
        {
            GSVector x   (particles[part][Bunch::x   ]);
            GSVector xp  (particles[part][Bunch::xp  ]);
            GSVector y   (particles[part][Bunch::y   ]);
            GSVector yp  (particles[part][Bunch::yp  ]);
            GSVector cdt (particles[part][Bunch::cdt ]);
            GSVector dpop(particles[part][Bunch::dpop]);

            // upstream slot
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

            // upstream edge
            FF_algorithm::edge_unit(y, yp, edge_k);

            // bend
            // FF_algorithm::dipole_unit(x, xp, y, yp, cdt, dpop, length, k[0]);
            FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop, 
                    dphi, strength, pref, m, reference_cdt, 
                    phase, term);

            // downstream edge
            FF_algorithm::edge_unit(y, yp, edge_k);

            // downstream slot
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

            particles[part][Bunch::x]   = x;
            particles[part][Bunch::xp]  = xp;
            particles[part][Bunch::y]   = y;
            particles[part][Bunch::yp]  = yp;
            particles[part][Bunch::cdt] = cdt;
 
        }

        for (int part = block_last; part < local_num; ++part) 
#else
        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part) 
#endif
        {
            double x   (particles[part][Bunch::x   ]);
            double xp  (particles[part][Bunch::xp  ]);
            double y   (particles[part][Bunch::y   ]);
            double yp  (particles[part][Bunch::yp  ]);
            double cdt (particles[part][Bunch::cdt ]);
            double dpop(particles[part][Bunch::dpop]);

            // upstream slot
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

            // upstream edge
            FF_algorithm::edge_unit(y, yp, edge_k);

            // bend
            // FF_algorithm::dipole_unit(x, xp, y, yp, cdt, dpop, length, k[0]);
            FF_algorithm::bend_unit(x, xp, y, yp, cdt, dpop, 
                    dphi, strength, pref, m, reference_cdt, 
                    phase, term);

            // downstream edge
            FF_algorithm::edge_unit(y, yp, edge_k);

            // downstream slot
            FF_algorithm::slot_unit(x, xp, y, yp, cdt, dpop, ct, st, pref, m);

            particles[part][Bunch::x]   = x;
            particles[part][Bunch::xp]  = xp;
            particles[part][Bunch::y]   = y;
            particles[part][Bunch::yp]  = yp;
            particles[part][Bunch::cdt] = cdt;
        }
    }
    else
    {
        double step_length = length/steps;
        double step_strength[6] = { k[0]*step_length, k[1]*step_length,
                                    k[2]*step_length, k[3]*step_length,
                                    k[4]*step_length, k[5]*step_length };

        double reference_cdt = get_reference_cdt(length, k, bunch.get_reference_particle());
        double step_reference_cdt = reference_cdt/steps;

        // with combined high order function, use yoshida approximation
        #pragma omp parallel for
        for (int part = 0; part < local_num; ++part) 
        {
            double x   (particles[part][Bunch::x   ]);
            double xp  (particles[part][Bunch::xp  ]);
            double y   (particles[part][Bunch::y   ]);
            double yp  (particles[part][Bunch::yp  ]);
            double cdt (particles[part][Bunch::cdt ]);
            double dpop(particles[part][Bunch::dpop]);

            FF_algorithm::yoshida<double, FF_algorithm::thin_rbend_unit<double>, 4, 3 >
                ( x, xp, y, yp, cdt, dpop,
                  reference_momentum, m,
                  step_reference_cdt,
                  step_length, step_strength, steps );

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
    FF_rbend::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_rbend::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_rbend::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_rbend::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_rbend::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_rbend::~FF_rbend()
{

}

