#include "ff_vkicker.h"
#include "ff_algorithm.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/utils/gsvector.h"

double FF_vkicker::get_reference_cdt(double length, double k, int steps, Reference_particle &reference_particle)
{
    double pref = reference_particle.get_momentum();
    double m = reference_particle.get_mass();
    double step_length = length / steps;
    // for 0 length, hk,vk is the total strength of the kick
    // for >0 length, hk,vk is the strength/length of the kick
    double step_strength[2] = { 0.0, k * step_length };

    double x(reference_particle.get_state()[Bunch::x]);
    double xp(reference_particle.get_state()[Bunch::xp]);
    double y(reference_particle.get_state()[Bunch::y]);
    double yp(reference_particle.get_state()[Bunch::yp]);
    double cdt(0.0);
    double dpop(reference_particle.get_state()[Bunch::dpop]);

#if 0
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
    FF_algorithm::thin_kicker_unit(yp, k);
    FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
#endif

    if ( close_to_zero(length) )
    {
        FF_algorithm::thin_kicker_unit(yp, k);
    }
    else
    {
        FF_algorithm::yoshida6<double, FF_algorithm::thin_kicker_unit<double>, 1>
            ( x, xp, y, yp, cdt, dpop,
              pref, m, 0.0,
              step_length, step_strength, steps );
    }

    reference_particle.set_state(x, xp, y, yp, cdt, dpop);

    return cdt;
}

void FF_vkicker::apply(Lattice_element_slice const& slice, JetParticle& jet_particle)
{
    throw std::runtime_error("libFF::vkicker for JetParticles has yet to be implemented");

#if 0
    double length = slice.get_right() - slice.get_left();

    double l = slice.get_lattice_element().get_double_attribute("l", 0.0);
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

    if ( close_to_zero(length) )
    {
        FF_algorithm::thin_kicker_unit(yp, k);
    }
    else
    {
        k = k * length / l;

        double pref = jet_particle.ReferenceMomentum();
        double m = jet_particle.Mass();

        Particle chef_particle(jet_particle);
        Reference_particle reference_particle(
                    chef_particle_to_reference_particle(chef_particle));

        double reference_cdt = get_reference_cdt(length, k, reference_particle);
        double step_reference_cdt = reference_cdt * 0.5;
        double step_length = length * 0.5;

        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
        FF_algorithm::thin_kicker_unit(yp, k);
        FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, 0.0);
    }
#endif
}

void FF_vkicker::apply(Lattice_element_slice const& slice, Bunch& bunch)
{
    double length = slice.get_right() - slice.get_left();

    double l = slice.get_lattice_element().get_double_attribute("l", 0.0);
    double k = slice.get_lattice_element().get_double_attribute("kick");

    Reference_particle       & ref_lattice = bunch.get_design_reference_particle();
    Reference_particle const & ref_bunch   = bunch.get_reference_particle();

    double plattice = ref_lattice.get_momentum();
    double pbunch = ref_bunch.get_momentum();
    // scale is to scale the kick strength defined relative to the lattice momentum to
    // the scale of the bunch particles defined relative to the bunch momentum
    double scale = plattice/pbunch;

    // kick strength is defined as momentum change/reference momentum
    k = k * ref_bunch.get_charge() / ref_lattice.get_charge();

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

    if ( close_to_zero(length) )
    {
        // the reference time is calculated with the design reference particle which
        // relative to plattice.
        // also update the reference particle
        double reference_cdt = get_reference_cdt(0.0, k, 1, ref_lattice);

        #pragma omp parallel for
        for (int part = 0; part < block_last; part += gsvsize)
        {
            GSVector yp(&ypa[part]);
            FF_algorithm::thin_kicker_unit(yp, k*scale);
            yp.store(&ypa[part]);
        }

        for (int part = block_last; part < local_num; ++part)
        {
            double yp(ypa[part]);
            FF_algorithm::thin_kicker_unit(yp, k*scale);
            ypa[part] = yp;
        }
    }
    else
    {
        // yoshida steps
        steps = (int)slice.get_lattice_element().get_double_attribute("yoshida_steps", 4.0);

        // strength per unit length
        k = k / l;

        double pref = bunch.get_reference_particle().get_momentum();
        double m = bunch.get_mass();
        double reference_cdt = get_reference_cdt(length, k, steps, ref_lattice);

        double step_reference_cdt = reference_cdt / steps;
        double step_length = length / steps;
        double step_strength[2] = { 0.0, k * step_length * scale};

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
            // simple drift - kick - drift scheme. this is how CHEF does it
            FF_algorithm::drift_unit(x, xp, y, yp, cdt, dpop, step_length, pref, m, step_reference_cdt);
            FF_algorithm::thin_kicker_unit(yp, k);
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
            FF_algorithm::thin_kicker_unit(yp, k);
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
    FF_vkicker::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(FF_element);
    }

template
void
FF_vkicker::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
FF_vkicker::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
FF_vkicker::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
FF_vkicker::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_vkicker::~FF_vkicker()
{

}

