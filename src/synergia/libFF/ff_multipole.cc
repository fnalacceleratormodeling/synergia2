#include "ff_multipole.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"

FF_multipole::FF_multipole()
{

}

void FF_multipole::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("FF_multipole::apply(JetParticle) not implemented");

#if 0
    double length = slice.get_right() - slice.get_left();
    double strength = 0.0;

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

    multipole_unit(x, xp, y, yp, cdt, dpop,
                   length, strength,
                   reference_momentum, m, reference_brho);
#endif
}

void FF_multipole::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    if (length > 0.0)
        throw std::runtime_error("FF_multipole::apply() cannot deal with thick elements");

    std::vector<double> const & knl = 
        slice.get_lattice_element().get_vector_attribute("knl");

    std::vector<double> const & ksl = 
        slice.get_lattice_element().get_vector_attribute("ksl");

    if (knl.size() != ksl.size())
        throw std::runtime_error("FF_multipole:apply() different array size of knl and ksl");

    double kL[2];

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    for (int part = 0; part < local_num; ++part) {
        double x   (particles[part][Bunch::x   ]);
        double xp  (particles[part][Bunch::xp  ]);
        double y   (particles[part][Bunch::y   ]);
        double yp  (particles[part][Bunch::yp  ]);

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

        // higher orders
        for (int n = 3; n < knl.size(); ++n) {
            if (knl[n] && ksl[n]) {
                kL[0] = knl[n]; kL[1] = ksl[n];
                FF_algorithm::thin_magnet_unit(x, xp, y, yp, kL, n);
            }
        }

        particles[part][Bunch::xp] = xp;
        particles[part][Bunch::yp] = yp;
    }
}

template<class Archive>
    void
    FF_multipole::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_multipole::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_multipole::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_multipole::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_multipole::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_multipole::~FF_multipole()
{

}

