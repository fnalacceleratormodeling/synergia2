#include "ff_kicker.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

double FF_kicker::get_reference_cdt(double length, double hk, double vk, Reference_particle &reference_particle) 
{
    double reference_cdt;

    double pref = reference_particle.get_momentum();
    double m = reference_particle.get_mass();
    double step_length = length / steps;
    double step_strength[2] = { hk * step_length, vk * step_length };

    double x   (reference_particle.get_state()[Bunch::x]);
    double xp  (reference_particle.get_state()[Bunch::xp]);
    double y   (reference_particle.get_state()[Bunch::y]);
    double yp  (reference_particle.get_state()[Bunch::yp]);
    double cdt (reference_particle.get_state()[Bunch::cdt]);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

    double cdt_orig = cdt;

#if 0
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
    FF_algorithm::thin_kicker_unit(xp, yp, hk, vk);
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
#endif

    if ( close_to_zero(length) ) 
    {
        FF_algorithm::thin_kicker_unit(xp, hk);
        FF_algorithm::thin_kicker_unit(yp, vk);
    }
    else
    {
        FF_algorithm::yoshida6<double, FF_algorithm::thin_kicker_unit<double>, 1>
            ( x, xp, y, yp, cdt, dpop,
              pref, m, 0.0, 
              step_length, step_strength, steps );
    }

    reference_particle.set_state(x, xp, y, yp, 0.0, dpop);
    reference_cdt = cdt - cdt_orig;

    return reference_cdt;
}

void FF_kicker::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF::kicker for JetParticles has yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

    double  l = slice.get_lattice_element().get_double_attribute("l", 0.0);
    double hk = slice.get_lattice_element().get_double_attribute("hkick");
    double vk = slice.get_lattice_element().get_double_attribute("vkick");

    typedef PropagatorTraits<JetParticle>::State_t State_t;
    typedef PropagatorTraits<JetParticle>::Component_t Component_t;

    State_t& state = jet_particle.State();

    Component_t & x(state[Chef::x]);
    Component_t & xp(state[Chef::xp]);
    Component_t & y(state[Chef::y]);
    Component_t & yp(state[Chef::yp]);
    Component_t & cdt(state[Chef::cdt]);
    Component_t & dpop(state[Chef::dpop]);

    if ( close_to_zero(l) ) 
    {
        double k[2] = { hk, vk };
        FF_algorithm::thin_kicker_unit(x, xp, y, yp, k);
    } 
    else 
    {
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        double k[2] = { hk*length/l, vk*length/l };

        double pref = jet_particle.ReferenceMomentum();
        double m = jet_particle.Mass();

        Particle chef_particle(jet_particle);
        Reference_particle reference_particle(
                    chef_particle_to_reference_particle(chef_particle));

        double reference_cdt = get_reference_cdt(length, hk, vk, reference_particle);
        double step_reference_cdt = reference_cdt / steps;
        double step_length = length / steps;

#if 0
        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
        FF_algorithm::thin_kicker_unit(xp, yp, hk, vk);
        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
#endif

        FF_algorithm::yoshida6<TJet<double>, FF_algorithm::thin_quadrupole_unit<TJet<double> >, 1 >
                ( x, xp, y, yp, cdt, dpop,
                  pref, m, step_reference_cdt,
                  step_length, k, steps );
 
    }
#endif
}

void FF_kicker::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    double  l = slice.get_lattice_element().get_double_attribute("l", 0.0);
    double hk = slice.get_lattice_element().get_double_attribute("hkick");
    double vk = slice.get_lattice_element().get_double_attribute("vkick");

    Reference_particle const & ref_lattice = get_ref_particle_from_slice(slice);
    Reference_particle const & ref_bunch   = bunch.get_reference_particle();

    hk = hk * ref_bunch.get_charge() / ref_lattice.get_charge();
    vk = vk * ref_bunch.get_charge() / ref_lattice.get_charge();

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    double * RESTRICT xa;
    double * RESTRICT xpa;
    double * RESTRICT ya;
    double * RESTRICT ypa;
    double * RESTRICT cdta;
    double * RESTRICT dpopa;

    bunch.set_arrays(xa, xpa, ya, ypa, cdta, dpopa);

    const int gsvsize = GSVector::size();
    const int num_blocks = local_num / gsvsize;
    const int block_last = num_blocks * gsvsize;

    if ( close_to_zero(l) ) 
    {
        GSVector gdummy(1.0);
        double   ddummy(1.0);

        double   k[2] = { hk, vk };

        // only to update the reference particle
        double reference_cdt = get_reference_cdt(0.0, hk, vk, bunch.get_reference_particle());

        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize) 
        {
            GSVector xp(&xpa[part]);
            GSVector yp(&ypa[part]);

            FF_algorithm::thin_kicker_unit(gdummy, xp, gdummy, yp, k);

            xp.store(&xpa[part]);
            yp.store(&ypa[part]);
        }

        for (int part = block_last; part < local_num; ++part) 
        {
            double xp(xpa[part]);
            double yp(ypa[part]);

            FF_algorithm::thin_kicker_unit(ddummy, xp, ddummy, yp, k);

            xpa[part] = xp;
            ypa[part] = yp;
        }
    } 
    else 
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // strength per unit length
        hk = hk / l;
        vk = vk / l;

        double pref = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, hk, vk, bunch.get_reference_particle());

        double step_reference_cdt = reference_cdt / steps;
        double step_length = length / steps;
        double step_strength[2] = { hk * step_length, vk * step_length };

        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize) 
        {
            GSVector x(&xa[part]);
            GSVector xp(&xpa[part]);
            GSVector y(&ya[part]);
            GSVector yp(&ypa[part]);
            GSVector cdt(&cdta[part]);
            GSVector dpop(&dpopa[part]);

#if 0
            // simple drift - kick - drift scheme, this is how CHEF does it
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
            FF_algorithm::thin_kicker_unit(xp, yp, hk, vk);
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
#endif

            FF_algorithm::yoshida6<GSVector, FF_algorithm::thin_kicker_unit<GSVector>, 1>
                ( x, xp, y, yp, cdt, dpop,
                  pref, m, step_reference_cdt, 
                  step_length, step_strength, steps );

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

#if 0
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
            FF_algorithm::thin_kicker_unit(xp, yp, hk, vk);
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
#endif

            FF_algorithm::yoshida6<double, FF_algorithm::thin_kicker_unit<double>, 1>
                ( x, xp, y, yp, cdt, dpop,
                  pref, m, step_reference_cdt, 
                  step_length, step_strength, steps );

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
    FF_kicker::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_kicker::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_kicker::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_kicker::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_kicker::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_kicker::~FF_kicker()
{

}

