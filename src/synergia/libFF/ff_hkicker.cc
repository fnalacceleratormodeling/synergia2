#include "ff_hkicker.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

double FF_hkicker::get_reference_cdt(double length, double k, Reference_particle &reference_particle) 
{
    if (length == 0) return 0.0;

    double reference_cdt;

    double pref = reference_particle.get_momentum();
    double m = reference_particle.get_mass();
    double step_length = length * 0.5;

    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(reference_particle.get_state()[Bunch::cdt]);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    double cdt_orig = cdt;

    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
    FF_algorithm::thin_hkicker_unit(xp, k);
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);

    reference_cdt = cdt - cdt_orig;

    return reference_cdt;
}

void FF_hkicker::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    double length = slice.get_right() - slice.get_left();

    double k = slice.get_lattice_element().get_double_attribute("kick");

    typedef PropagatorTraits<JetParticle>::State_t State_t;
    typedef PropagatorTraits<JetParticle>::Component_t Component_t;

    State_t& state = jet_particle.State();

    Component_t & x(state[Chef::x]);
    Component_t & xp(state[Chef::xp]);
    Component_t & y(state[Chef::y]);
    Component_t & yp(state[Chef::yp]);
    Component_t & cdt(state[Chef::cdt]);
    Component_t & dpop(state[Chef::dpop]);

    double pref = jet_particle.ReferenceMomentum();
    double m = jet_particle.Mass();

    Particle chef_particle(jet_particle);
    Reference_particle reference_particle(
                chef_particle_to_reference_particle(chef_particle));

    double reference_cdt = get_reference_cdt(length, k, reference_particle);
    double step_reference_cdt = reference_cdt * 0.5;
    double step_length = length * 0.5;

    if (length == 0.0) 
    {
        FF_algorithm::thin_hkicker_unit(yp, k);
    } 
    else 
    {
        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
        FF_algorithm::thin_hkicker_unit(xp, k);
        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
    }
}

void FF_hkicker::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();
    double k = slice.get_lattice_element().get_double_attribute("kick");

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double * RESTRICT xa;
    double * RESTRICT xpa;
    double * RESTRICT ya;
    double * RESTRICT ypa;
    double * RESTRICT cdta;
    double * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int num_blocks = local_num / GSVector::size;
    const int block_last = num_blocks * GSVector::size;

    if (length == 0.0) 
    {
        for (int part = 0; part < block_last; part += GSVector::size) 
        {
            GSVector xp(&xpa[part]);
            FF_algorithm::thin_hkicker_unit(xp, k);
            xp.store(&xpa[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double xp(xpa[part]);
            FF_algorithm::thin_hkicker_unit(xp, k);
            xpa[part] = xp;
        }
    } 
    else 
    {
        double pref = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k, bunch.get_reference_particle());

        double step_reference_cdt = reference_cdt * 0.5;
        double step_length = length * 0.5;

        for (int part = 0; part < block_last; part += GSVector::size) 
        {
            GSVector x(&xa[part]);
            GSVector xp(&xpa[part]);
            GSVector y(&ya[part]);
            GSVector yp(&ypa[part]);
            GSVector cdt(&cdta[part]);
            GSVector dpop(&dpopa[part]);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
            FF_algorithm::thin_hkicker_unit(xp, k);
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);

            x.store(&xa[part]);
            xp.store(&xpa[part]);
            y.store(&ya[part]);
            yp.store(&ypa[part]);
            cdt.store(&cdta[part]);
            dpop.store(&dpopa[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double x(xa[part]);
            double xp(xpa[part]);
            double y(ya[part]);
            double yp(ypa[part]);
            double cdt(cdta[part]);
            double dpop(dpopa[part]);

            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
            FF_algorithm::thin_hkicker_unit(xp, k);
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);

            xa[part] = x;
            xpa[part] = xp;
            ya[part] = y;
            ypa[part] = yp;
            cdta[part] = cdt;
            dpopa[part] = dpop;
        }

        bunch.get_reference_particle().increment_trajectory(length);
    }
}

template<class Archive>
    void
    FF_hkicker::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_hkicker::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_hkicker::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_hkicker::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_hkicker::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_hkicker::~FF_hkicker()
{

}

