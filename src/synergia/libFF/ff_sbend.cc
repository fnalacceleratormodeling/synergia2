#include "ff_sbend.h"
#include "synergia/lattice/chef_utils.h"

FF_sbend::FF_sbend()
{

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
}

void FF_sbend::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();
    double angle = slice.get_lattice_element().get_double_attribute("angle");
    double l     = slice.get_lattice_element().get_double_attribute("l");

    double cos_angle = cos(angle);
    double sin_angle = sin(angle);

    double e1 = 0.0;
    double e2 = 0.0;

    if (slice.get_lattice_element().has_double_attribute("e1"))
        e1 = slice.get_lattice_element().get_double_attribute("e1");
    if (slice.get_lattice_element().has_double_attribute("e2"))
        e2 = slice.get_lattice_element().get_double_attribute("e2");


    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double reference_momentum = bunch.get_reference_particle().get_momentum();
    double reference_brho     = reference_momentum / PH_CNV_brho_to_p;
    double m = bunch.get_mass();

    double strength = reference_brho * angle / l;

    double reference_cdt = get_reference_cdt(length, angle, strength,
                                             bunch.get_reference_particle());

    double psi = angle - (e1 + e2);
    std::complex<double> phase = std::exp( std::complex<double>(0.0, psi) );
    std::complex<double> term = std::complex<double>(0.0, length / angle) *
                                std::complex<double>(1.0 - cos_angle, - sin_angle) *
                                std::complex<double>(cos(e2), -sin(e2));


    double pref = reference_momentum;

    double ce1 = cos(-e1);
    double se1 = sin(-e1);
    double ce2 = cos(-e2);
    double se2 = sin(-e2);


    for (int part = 0; part < local_num; ++part) {
        double x   (particles[part][Bunch::x   ]);
        double xp  (particles[part][Bunch::xp  ]);
        double y   (particles[part][Bunch::y   ]);
        double yp  (particles[part][Bunch::yp  ]);
        double cdt (particles[part][Bunch::cdt ]);
        double dpop(particles[part][Bunch::dpop]);

            double r0 = x;
            double r1 = y;
            double r2 = 0.0;

            double zp = sqrt((dpop+1)*(dpop+1) - xp * xp - yp * yp);

            double p = (dpop+1) * pref;
            double e = sqrt(p*p + m*m);

            double b0 = xp * pref / e;
            double b1 = yp * pref / e;
            double b2 = zp * pref / e;

            double bp = se1 * b0 + ce1 * b2;

            double tau = -x * se1 / bp;

            std::cout << std::setprecision(10);

            std::cout << "state = " << x << ", " << y << ", " << xp << ", " << yp <<"\n";

            std::cout << "myb = " << b0 << ", " << b1 << ", " << b2 << "\n";
            std::cout << "mybp = " << bp << "\n";
            std::cout << "mytau = " << tau << "\n";
            std::cout << "cos(t) = " << ce1 << ", sin(t) = " << se2 << "\n";
            std::cout << "p = " << xp << ", " << yp << ", " << zp << "\n";

            r0 = x + tau * b0;
            r1 = y + tau * b1;
            r2 = 0 + tau * b2;

            x = r0 * ce1 - r2 * se1;
            y = r1;

            xp = xp * ce1 - zp * se1;
            yp = yp;

            std::cout << "state after slot 1 = " << x << ", " << y << ", " << xp << ", " << yp <<"\n";

        sbend_unit2(x, xp, y, yp, cdt, dpop,
                   length, angle, strength,
                   reference_momentum, m, reference_cdt, phase, term);

            std::cout << "state after sbend  = " << x << ", " << y << ", " << xp << ", " << yp <<"\n";

            b0 = xp * pref / e;
            b1 = yp * pref / e;

            bp = se2 * b0 + ce2 * b2;
            tau = -x * se2 / bp;

            std::cout << "myb = " << b0 << ", " << b1 << ", " << b2 << "\n";
            std::cout << "mybp = " << bp << "\n";
            std::cout << "mytau = " << tau << "\n";
            std::cout << "cos(t) = " << ce2 << ", sin(t) = " << se2 << "\n";
            std::cout << "p = " << xp << ", " << yp << ", " << zp << "\n";

            r0 = x + tau * b0;
            r1 = y + tau * b1;
            r2 = 0 + tau * b2;

            x = r0 * ce2 - r2 * se2;
            y = r1;

            xp = xp * ce2 - zp * se2;

            std::cout << "px = " << xp << ", py = " << yp << "\n";


        particles[part][Bunch::x]  = x;
        particles[part][Bunch::xp] = xp;
        particles[part][Bunch::y]  = y;
        particles[part][Bunch::cdt] = cdt;
    }

    bunch.get_reference_particle().increment_trajectory(length);

#if 0
    double length = slice.get_right() - slice.get_left();
    double angle = slice.get_lattice_element().get_double_attribute("angle");

    double cos_angle = cos(angle);
    double sin_angle = sin(angle);

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double reference_momentum = bunch.get_reference_particle().get_momentum();
    double reference_brho     = reference_momentum / PH_CNV_brho_to_p;
    double m = bunch.get_mass();

    std::cout << "brho = " << reference_brho << "\n";

    for (int part = 0; part < local_num; ++part) {
        double x   (particles[part][Bunch::x   ]);
        double xp  (particles[part][Bunch::xp  ]);
        double y   (particles[part][Bunch::y   ]);
        double yp  (particles[part][Bunch::yp  ]);
        double cdt (particles[part][Bunch::cdt ]);
        double dpop(particles[part][Bunch::dpop]);

        sbend_unit(x, xp, y, yp, cdt, dpop,
                   length, cos_angle, sin_angle,
                   reference_momentum, m, reference_brho);

        particles[part][Bunch::x]  = x;
        particles[part][Bunch::xp] = xp;
        particles[part][Bunch::y]  = y;
        particles[part][Bunch::cdt] = cdt;
    }

    bunch.get_reference_particle().increment_trajectory(length);
#endif
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

