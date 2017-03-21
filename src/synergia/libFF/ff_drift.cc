#include "ff_drift.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

FF_drift::FF_drift()
{
}

double FF_drift::get_reference_cdt(double length, Reference_particle const& reference_particle)
{
    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(reference_particle.get_state()[Bunch::cdt]);
    double dpop(reference_particle.get_state()[Bunch::dpop]);
    double reference_momentum = reference_particle.get_momentum();
    double m = reference_particle.get_mass();

    double cdt_orig = cdt;
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m, 0.0);

    return cdt - cdt_orig;
}

void FF_drift::apply(Lattice_element_slice const& slice,
                     Trigon_particle_t & trigon_particle)
{
    const double  length = slice.get_right() - slice.get_left();
    const double    mass = trigon_particle.get_reference_particle().get_mass();
    Reference_particle const& ref_b = trigon_particle.get_reference_particle();
    const double ref_cdt = get_reference_cdt(length, ref_b);


    const double   ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    Trigon_particle_t::Component_t& x(trigon_particle.get_state()[Bunch::x]);
    Trigon_particle_t::Component_t& xp(trigon_particle.get_state()[Bunch::xp]);
    Trigon_particle_t::Component_t& y(trigon_particle.get_state()[Bunch::y]);
    Trigon_particle_t::Component_t& yp(trigon_particle.get_state()[Bunch::yp]);
    Trigon_particle_t::Component_t& cdt(trigon_particle.get_state()[Bunch::cdt]);
    Trigon_particle_t::Component_t& dpop(trigon_particle.get_state()[Bunch::dpop]);

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);
}

void FF_drift::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    const double  length = slice.get_right() - slice.get_left();
    const double    mass = jet_particle.Mass();
    Particle particle(jet_particle);
    Reference_particle ref_b(chef_particle_to_reference_particle(particle));
    const double ref_cdt = get_reference_cdt(length, ref_b);


    const double   ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    Jet& x(jet_particle.State()[Particle::xIndex]);
    Jet& xp(jet_particle.State()[Particle::npxIndex]);
    Jet& y(jet_particle.State()[Particle::yIndex]);
    Jet& yp(jet_particle.State()[Particle::npyIndex]);
    Jet& cdt(jet_particle.State()[Particle::cdtIndex]);
    Jet& dpop(jet_particle.State()[Particle::ndpIndex]);

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);
}

void FF_drift::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const int  local_num = bunch.get_local_num();
    const double    mass = bunch.get_mass();
    const double ref_cdt = get_reference_cdt(length, bunch.get_reference_particle());

    Reference_particle & ref_b = bunch.get_reference_particle();
    const double   ref_p = ref_b.get_momentum() * (1.0 + ref_b.get_state()[Bunch::dpop]);

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

    #pragma omp parallel for
    for (int part = 0; part < block_last; part += gsvsize) 
    {
        GSVector x(&xa[part]);
        GSVector xp(&xpa[part]);
        GSVector y(&ya[part]);
        GSVector yp(&ypa[part]);
        GSVector cdt(&cdta[part]);
        GSVector dpop(&dpopa[part]);

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

        x.store(&xa[part]);
        y.store(&ya[part]);
        cdt.store(&cdta[part]);
    }

    for (int part = block_last; part < local_num; ++part) 
    {
        double x(xa[part]);
        double xp(xpa[part]);
        double y(ya[part]);
        double yp(ypa[part]);
        double cdt(cdta[part]);
        double dpop(dpopa[part]);

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, ref_p, mass, ref_cdt);

        xa[part] = x;
        ya[part] = y;
        cdta[part] = cdt;
    }

    bunch.get_reference_particle().increment_trajectory(length);
}

template<class Archive>
    void
    FF_drift::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_drift::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_drift::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_drift::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_drift::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_drift::~FF_drift()
{

}
