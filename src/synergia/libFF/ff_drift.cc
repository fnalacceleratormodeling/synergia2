#include "ff_drift.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"
#include "synergia/utils/logger.h"

FF_drift::FF_drift()
{
}

double FF_drift::get_reference_cdt(double length, Reference_particle & reference_particle)
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

void FF_drift::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    double length = slice.get_right() - slice.get_left();

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
    double m = jet_particle.Mass();

    Particle chef_particle(jet_particle);
    Reference_particle reference_particle(
                chef_particle_to_reference_particle(chef_particle));
    double reference_cdt = get_reference_cdt(length, reference_particle);

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, length, reference_momentum, m,
               reference_cdt);
}

void FF_drift::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    const double  length = slice.get_right() - slice.get_left();
    const int  local_num = bunch.get_local_num();
    const double   ref_p = bunch.get_reference_particle().get_momentum();
    const double    mass = bunch.get_mass();
    const double ref_cdt = get_reference_cdt(length, bunch.get_reference_particle());

    double * RESTRICT xa, * RESTRICT xpa;
    double * RESTRICT ya, * RESTRICT ypa;
    double * RESTRICT cdta, * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int num_blocks = local_num / GSVector::size;
    const int block_last = num_blocks * GSVector::size;

    #pragma omp parallel for
    for (int part = 0; part < block_last; part += GSVector::size) 
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
