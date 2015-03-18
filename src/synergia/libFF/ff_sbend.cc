#include "ff_sbend.h"
#include "synergia/lattice/chef_utils.h"

FF_sbend::FF_sbend()
{

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

    double cos_angle = cos(angle);
    double sin_angle = sin(angle);

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double reference_momentum = bunch.get_reference_particle().get_momentum();
    double reference_brho     = reference_momentum / PH_CNV_brho_to_p;
    double m = bunch.get_mass();

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

